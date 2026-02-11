# LEASH Comparisons

Canonical comparison document for LEASH.

## Scope

This file centralizes comparisons across:

- LEASH topologies (width and merge shape).
- Data movement models (holepunch vs S3-mediated).
- S3 reader strategies.
- Reliability and performance tradeoffs.

For implementation details, use:

- `LEASH.md`
- `S3_DIRECT_STREAMING_METHODS.md`

## 1) Topology Comparison: Width 2 vs Width 4

| Dimension | Width 2 | Width 4 |
| --- | --- | --- |
| Split fanout | 2 workers | 4 workers |
| Merge shape | Single merge stage | Two-level merge tree |
| Intermediate EC2 merge nodes | 0 | 2 |
| Coordination complexity | Lower | Higher |
| Typical use | Simpler/medium workloads | Higher parallelism workloads |

Key insight:

- Width 4 increases parallel work but adds merge stages and coordination overhead.

## 2) Data Path Comparison

| Model | Path | Main advantage | Main risk |
| --- | --- | --- | --- |
| Classic LEASH transport | `S3 -> EC2 split -> pashlib -> Lambda -> pashlib -> EC2 merge -> S3` | Mature pipeline shape | Holepunch timing/coordination sensitivity |
| S3 direct readers | `S3 -> Lambda readers -> EC2 merge -> S3` | Removes EC2 full-file download/split | Boundary correctness and stream backpressure tuning |
| Manual S3 orchestration | `S3 <-> Lambda via explicit orchestrator scripts` | Operational reliability and debuggability | Not the full automatic LEASH transport path |

## 3) S3 Reader Strategy Comparison

| Strategy ID | Script | Accuracy model | Coordination | Recommended use |
| --- | --- | --- | --- | --- |
| `s3_approx_correction` | `aws/s3-chunk-reader-approx-correction.py` | Lambda-side boundary correction | Independent workers | Default choice |
| `s3_smart_prealigned` | `aws/s3-chunk-reader-smart-prealigned.py` | EC2 prealigned ranges | Independent workers | Correctness validation/debug |
| `s3_approx_tail_coord` | `aws/s3-chunk-reader-approx-tail-coord.py` | Tail correction via predecessor info | Chain dependencies | Research/legacy only |

Decision rule:

1. Start with `s3_approx_correction`.
2. If investigating boundary bugs, compare with `s3_smart_prealigned`.
3. Use `s3_approx_tail_coord` only for explicit coordination experiments.

## 4) Reliability Comparison

Historical operational behavior:

| Approach | Typical behavior |
| --- | --- |
| Holepunch transport path | Sensitive to startup timing and NAT traversal conditions |
| S3-mediated manual orchestration | More predictable due to explicit S3 artifacts and synchronous worker flow |

Interpretation:

- If the goal is stable debugging and deterministic repro, S3-mediated orchestration can be easier to operate.
- If the goal is end-to-end LEASH transport validation, holepunch path is still the primary architecture to test.

## 5) Historical Performance Snapshot

These values are historical and workload-dependent. Re-run benchmarks for current numbers.

### S3 strategy snapshot

| Workload | Baseline EC2 split | Smart prealigned | Tail coordination |
| --- | --- | --- | --- |
| 1MB text | 4-7s | 2.5-8s | ~30s |
| 1GB text | 20-47s | ~30s | ~28s |
| 20GB text | 3-4m | 5.5-6m | frequent timeouts |

### Manual S3 orchestrator snapshot (1MB, historical)

| Variant | Total time |
| --- | --- |
| Full download/split | ~2.11s |
| No-split | ~1.71s |
| Byte-range | ~1.78s |

## 6) How To Reproduce Comparisons

Benchmark matrix (recommended):

```bash
cd /home/ubuntu/pash/evaluation/benchmarks
./run-leash-benchmark-matrix.sh --help
```

Manual S3 orchestration quick run:

```bash
cd /home/ubuntu/pash
python3 manual_s3_orchestrator_byte_ranges.py \
  --bucket <bucket> \
  --input oneliners/inputs/1M.txt \
  --output oneliners/outputs/result.txt \
  --workers 2
```

## 7) Source References

- Strategy selection and wiring: `compiler/serverless/ir_helper.py`
- Reader strategy constants/mapping: `compiler/definitions/ir/nodes/serverless_remote_pipe.py`
- LEASH orchestration manager: `compiler/serverless/serverless_executor.py`
- Manual S3 orchestrators: `manual_s3_orchestrator.py`, `manual_s3_orchestrator_no_split.py`, `manual_s3_orchestrator_byte_ranges.py`
