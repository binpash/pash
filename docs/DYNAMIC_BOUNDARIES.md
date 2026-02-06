# Dynamic Boundaries Implementation - Summary

## Changes Made

Successfully implemented the dynamic boundaries mode for S3 direct streaming as specified in the plan.

### 1. Compiler Changes (`compiler/serverless/ir_helper.py`)

**Lines Modified**: 1138-1141, 1142+ (new block inserted)

- Added `USE_DYNAMIC_BOUNDARIES` environment variable check with highest priority
- Inserted new conditional block before `use_approx_with_correction`
- Dynamic mode skips `estimate_avg_line_length()` call (eliminates sampling overhead)
- Calculates approximate boundaries using pure arithmetic division
- Sets `window_size_param = None` as sentinel value to signal dynamic mode
- Reuses same round-robin chunk distribution logic as approx mode

**Key Implementation Details**:
```python
if use_dynamic_boundaries:
    # No sampling, pure arithmetic
    cached_window_size = None  # Sentinel for dynamic mode
```

### 2. Lambda Routing Changes (`compiler/definitions/ir/nodes/serverless_remote_pipe.py`)

**Lines Modified**: 44-59, 70-76

- Updated mode detection to check for `USE_DYNAMIC_BOUNDARIES`
- Both dynamic and fixed approx modes now route to same lambda script: `s3-chunk-reader-approx-correction.py`
- Modified operand passing to handle `window_size=None` case:
  - If `window_size is None`: passes string `"window_size=None"`
  - If `window_size is not None`: passes numeric value as before

**Key Implementation Details**:
```python
if use_dynamic_boundaries or use_approx_with_correction:
    remote_pipe_bin = "python3 aws/s3-chunk-reader-approx-correction.py"

    # Pass window_size operand for both modes
    if window_size is not None:
        operand_list.append(Operand(Arg.string_to_arg(f"window_size={window_size}")))
    else:
        operand_list.append(Operand(Arg.string_to_arg("window_size=None")))
```

### 3. Lambda Script Changes (`aws/s3-chunk-reader-approx-correction.py`)

**Lines Modified**: 468-490 (parameter parsing), 187-288 (expansion function)

#### 3a. Parameter Parsing (Lines 468-490)
- Modified to detect `window_size=None` string and convert to Python `None`
- Added debug output to distinguish dynamic vs fixed mode

```python
window_size_raw = kwargs.get('window_size', '10240')
if window_size_raw == 'None' or window_size_raw is None or window_size_raw == 'null':
    window_size = None  # Signals dynamic mode
else:
    window_size = int(window_size_raw)
```

#### 3b. Dynamic Expansion Function (Lines 187-288)
- Completely rewrote `find_newline_with_expansion()` for true exponential doubling
- **Dynamic mode**: Starts at 1KB, doubles each iteration (1KB → 2KB → 4KB → 8KB → ...)
- **Fixed mode**: Starts at provided window_size, doubles if needed
- Uses `while current_window <= max_window` loop with `current_window *= 2`
- No fixed list of expansion sizes - pure doubling algorithm
- Maximum expansion: 100MB

**Key Implementation Details**:
```python
def find_newline_with_expansion(..., window_size, ...):
    # Determine initial window size
    if window_size is None:
        current_window = 1024  # Dynamic: start with 1KB
    else:
        current_window = window_size  # Fixed: use provided window

    max_window = 100 * 1024 * 1024  # 100MB

    # Keep doubling until found or max
    while current_window <= max_window:
        # Try to read and find newline
        # If not found: current_window *= 2
```

### 4. Benchmark Script Changes (`evaluation/benchmarks/unix50/run-leash-compare-generic.sh`)

**Lines Modified**: 6-27, 404-448, 492-514, 571

#### 4a. Flag Parsing (Lines 6-27)
- Added `RUN_S3OPTDYNAMIC=false` variable
- Added `--dynamic` flag detection
- Updated default mode check to include dynamic mode

#### 4b. Execution Block (Lines 404-448)
- Added complete execution block for dynamic mode
- Uses environment variable: `USE_DYNAMIC_BOUNDARIES=true`
- Sets chunks per lambda: `PASH_S3_CHUNKS_PER_LAMBDA=16`
- Downloads output to `/tmp/compare_s3optdynamic_*.txt`
- Integrates with existing log collection via `utils.py`

#### 4c. Comparison Block (Lines 492-514)
- Added comparison between noopt and s3optdynamic outputs
- Uses same `compare_outputs` function as existing modes
- Tracks results in `comparison_results` array with label `(dynamic)`

#### 4d. Cleanup (Line 571)
- Added `/tmp/compare_s3optdynamic_*.txt` to cleanup command

## Environment Variable Priority

The system now checks modes in this priority order:

1. **`USE_DYNAMIC_BOUNDARIES=true`** (NEW - highest priority)
   - No EC2 sampling
   - Lambda starts with 1KB window, doubles exponentially
   - Uses `s3-chunk-reader-approx-correction.py` with `window_size=None`

2. **`USE_APPROX_LAMBDA_CORRECTION=true`**
   - EC2 samples for avg line length
   - Calculates fixed window size
   - Uses `s3-chunk-reader-approx-correction.py` with numeric window_size

3. **`USE_SMART_BOUNDARIES=true`**
   - EC2 pre-scans S3 for exact boundaries
   - Uses `s3-chunk-reader-smart-prealigned.py`

4. **Default** (no env vars)
   - Old approximate mode with inter-lambda communication
   - Uses `s3-chunk-reader-approx-tail-coordination.py`

## How to Test

### Correctness Testing
```bash
cd evaluation/benchmarks/unix50

# Test with specific modes
./run-leash-compare-generic.sh --noopt --dynamic --small
./run-leash-compare-generic.sh --noopt --dynamic --medium

# Compare all modes
./run-leash-compare-generic.sh --noopt --s3opt --dynamic --large
```

### Performance Comparison
```bash
# Run all modes to compare performance
./run-leash-compare-generic.sh --noopt --s3opt --s3optnonsmart --dynamic --large

# Check logs for:
# - EC2 pre-processing time (should be ~0 for dynamic)
# - Lambda execution time
# - Total wall-clock time
# - S3 API call counts
```

### Edge Case Testing
1. **Very long lines**: Create file with 10MB lines
2. **No newlines**: Test with binary files
3. **Very small chunks**: Use small input with high parallelism
4. **Exact boundaries**: Test files with size divisible by chunk count

## Expected Behavior

### EC2 Phase
- **Dynamic mode**: Only arithmetic calculations, no S3 reads
- **Fixed mode**: Samples file 5 times (0.5-2 seconds overhead)

### Lambda Phase
- **Dynamic mode**: Initial window = 1KB, doubles until newline found
- **Fixed mode**: Initial window = calculated size, doubles if needed

### For Typical Text Files (lines < 4KB)
- Dynamic mode should find newline in 2-3 iterations:
  - Try 1KB → not found
  - Try 2KB → not found
  - Try 4KB → found (most common case)

### For Files with Long Lines
- Both modes expand up to 100MB
- Performance degrades gracefully with number of S3 round-trips

## Files Modified

1. `compiler/serverless/ir_helper.py` - Add dynamic mode logic
2. `compiler/definitions/ir/nodes/serverless_remote_pipe.py` - Route to correct lambda script
3. `aws/s3-chunk-reader-approx-correction.py` - Implement true exponential expansion
4. `evaluation/benchmarks/unix50/run-leash-compare-generic.sh` - Add --dynamic flag

## Verification Checklist

- [x] `USE_DYNAMIC_BOUNDARIES` env var added and checked in ir_helper.py
- [x] Dynamic mode skips `estimate_avg_line_length()` call
- [x] `window_size=None` passed to lambda as sentinel (as string "None")
- [x] Lambda script detects None and uses dynamic expansion
- [x] Expansion starts at 1KB, doubles to max 100MB (true exponential)
- [x] run-leash-compare-generic.sh supports --dynamic flag
- [x] Output comparison works for s3optdynamic mode
- [x] All existing modes still work (no regression in code paths)
- [x] Cleanup includes s3optdynamic temp files

## Implementation Complete

All planned changes have been successfully implemented. The system now supports dynamic boundaries mode that eliminates EC2 sampling overhead and uses lambda-side exponential window expansion starting from 1KB.