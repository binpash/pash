# S3 Direct Streaming Methods

Canonical reference for S3-based reader strategies in LEASH.

## Scope

This document consolidates:

- S3 reader strategy evolution.
- Current strategy-to-script mapping used by the compiler/runtime.
- Manual S3 orchestration variants.
- Tuning knobs, known risks, and troubleshooting.

If you are deciding which S3 approach to use, start with the strategy table below.

## Current Strategy Map

The runtime supports three reader strategies. These names are the canonical strategy IDs used in code.

| Strategy ID | Reader script | Boundary handling | Coordination model | Status |
| --- | --- | --- | --- | --- |
| `s3_approx_correction` | `aws/s3-chunk-reader-approx-correction.py` | Approximate boundaries from EC2 + Lambda-side correction | Independent Lambdas | Preferred |
| `s3_smart_prealigned` | `aws/s3-chunk-reader-smart-prealigned.py` | EC2 pre-computes line-aligned ranges | Independent Lambdas | Supported |
| `s3_approx_tail_coord` | `aws/s3-chunk-reader-approx-tail-coord.py` | Approximate ranges corrected via predecessor coordination | Lambda-to-Lambda chain | Legacy/experimental |

Code references:

- Strategy constants and script mapping: `compiler/definitions/ir/nodes/serverless_remote_pipe.py`
- Per-lambda strategy selection: `compiler/serverless/ir_helper.py`

## How Strategy Selection Works

There are two paths:

1. Preferred path: `ir_helper.py` computes and passes `s3_reader_strategy` explicitly.
2. Backward-compatible path: `serverless_remote_pipe.py` derives strategy from environment flags.

Environment fallback behavior:

- Any of these flags set to `true` selects `s3_approx_correction`:
  - `USE_ADAPTIVE_BOUNDARIES`
  - `USE_DYNAMIC_BOUNDARIES`
  - `USE_ADAPTIVE_SIMPLE`
  - `USE_APPROX_LAMBDA_CORRECTION`
  - `USE_SINGLE_SHOT`
- If correction flags are off and `USE_SMART_BOUNDARIES=true`, select `s3_smart_prealigned`.
- Otherwise fallback to `s3_approx_tail_coord`.

## Method Evolution (Historical)

### Method 1: Basic S3 direct reads

- Idea: remove EC2 splitter overhead by letting Lambdas read directly.
- Result: correctness and performance issues.
- Status: superseded.

### Method 2: dgsh-tee eager buffering

- Idea: push more buffering into FIFOs/eager nodes.
- Result: memory pressure and sequential FIFO consumption effects.
- Status: superseded.

### Method 3: Tail coordination (DFR-style)

- Idea: each Lambda communicates boundary information to the next Lambda.
- Result: sequential dependency chain, high coordination overhead, timeout risk at scale.
- Status: kept as legacy strategy `s3_approx_tail_coord`.

### Method 4: Smart prealigned boundaries

- Idea: EC2 computes exact line-aligned boundaries before worker execution.
- Result: robust correctness, but EC2 pre-scan overhead can dominate for some workloads.
- Status: available as `s3_smart_prealigned`.

### Method 5: Approximate boundaries + Lambda correction

- Idea: EC2 does arithmetic-only chunking; each Lambda locally corrects start/end boundaries.
- Result: best tradeoff in most runs (parallelism preserved, low EC2 overhead, correct output).
- Status: preferred strategy `s3_approx_correction`.

## Manual S3 Orchestration Variants

These bypass holepunch transport and use S3 as the exchange medium.

| Variant | Script | Best use |
| --- | --- | --- |
| Full download/split/upload | `manual_s3_orchestrator.py` | Simple end-to-end baseline |
| No-split (pre-existing chunks) | `manual_s3_orchestrator_no_split.py` | Reusing already chunked inputs |
| Byte-range direct read | `manual_s3_orchestrator_byte_ranges.py` | Most efficient manual S3 path |

Related Lambda workers:

- `lambda_sort_worker.py`
- `lambda_sort_worker_byte_ranges.py`

Deployment helpers:

- `deploy_lambda_sort.sh`
- `deploy_lambda_sort_byte_ranges.sh`

## Tuning Knobs

### Core mode flags

- `USE_APPROX_LAMBDA_CORRECTION=true`
- `USE_SMART_BOUNDARIES=true`
- `USE_DYNAMIC_BOUNDARIES=true`
- `USE_ADAPTIVE_BOUNDARIES=true`
- `USE_ADAPTIVE_SIMPLE=true`
- `USE_SINGLE_SHOT=true`

### Correction/adaptive parameters

- `PASH_S3_CHUNKS_PER_LAMBDA`
- `PASH_ADAPTIVE_TARGET_LINES`
- `PASH_ADAPTIVE_RETRY_RISK`
- `PASH_ADAPTIVE_SAMPLE_KB`
- `PASH_ADAPTIVE_SAFETY_FACTOR`

Notes:

- Correction-based modes pass `window_size` operands to the reader.
- In dynamic mode, `window_size=None` is passed and Lambda expands probe windows progressively.
- Header handling is passed explicitly to readers (`write_headers=true|false`).

## Historical Performance Snapshot

These numbers are from earlier experiments and should be treated as directional only.

| Workload | Baseline EC2 split | Smart prealigned | Tail coordination |
| --- | --- | --- | --- |
| 1MB text | 4-7s | 2.5-8s | ~30s |
| 1GB text | 20-47s | ~30s | ~28s |
| 20GB text | 3-4m | 5.5-6m | frequent timeouts |

Additional historical result from method writeups:

- Approximate + correction reached around ~3m on 20GB in some runs.

Use the benchmark matrix runner for current numbers:

- `evaluation/benchmarks/run-leash-benchmark-matrix.sh`

## Known Failure Modes

### FIFO backpressure + S3 stream timeout

Symptom:

- `IncompleteReadError` while streaming from S3 into FIFO.

Typical cause:

- Downstream consumer is slower than producer; FIFO blocks long enough for S3 read timeout.

Mitigations:

- Reduce per-reader pressure via chunk sizing.
- Prefer correction strategy over tail coordination for less coordination blocking.
- Validate Lambda memory/ephemeral storage sizing for sorting workloads.

### Boundary mistakes (missing/duplicated lines)

Symptom:

- Output mismatch, dropped bytes, or duplicated boundary lines.

Mitigations:

- Ensure correction strategy performs previous-byte and newline alignment checks.
- Re-run with smart prealigned mode for strict boundary validation.

### Large-file slowdowns despite direct S3

Symptom:

- 20GB+ runs slower than EC2 split baseline.

Common hypotheses:

- S3 parallel read contention.
- Lambda memory/IO bottlenecks.
- Merge-phase pressure on coordinator.

## Recommended Defaults

- Start with `s3_approx_correction`.
- Use `s3_smart_prealigned` when validating strict boundary behavior or debugging correctness.
- Use `s3_approx_tail_coord` only when explicitly testing coordination-based behavior.
- If holepunch transport is unstable for a workflow, use manual S3 orchestration scripts.

## Related Docs

- LEASH architecture: `LEASH.md`
- Comparison hub: `LEASH_ANALYSIS.md`
- Refactor history: `S3_READER_REFACTOR_CHANGELOG.md`
