"""
S3 Boundary Detection - Single-Shot Mode

A stripped-down variant of adaptive_simple that issues exactly **one** S3 GET
for boundary setup, regardless of file size or shard count.  The single sample
is taken from the **middle** of the file (filesize // 2) rather than byte 0,
giving a more representative estimate than the shared helper would produce with
num_samples=1.  A higher safety factor (default 2.0 vs adaptive_simple's 1.5)
compensates for the reduced statistical coverage.

Environment Variables:
    - USE_SINGLE_SHOT=true: Enable this mode
    - PASH_SINGLE_SHOT_SAMPLE_KB: Sample size in KB (default 256)
    - PASH_SINGLE_SHOT_SAFETY_FACTOR: Multiplier on avg line length (default 2.0)

Functions:
    - calculate_single_shot_boundaries: 1 sample → window → arithmetic shard_ranges

Usage:
    from compiler.serverless.s3_single_shot import calculate_single_shot_boundaries

    shard_ranges, window_size = calculate_single_shot_boundaries(
        bucket="mybucket", key="file.txt", filesize=1000000,
        total_chunks=256, config=BoundaryConfig(), debug=True
    )  # window_size is an int (bytes)
"""

import boto3, sys, os
from typing import List, Tuple

sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
from serverless.s3_config import S3BoundaryConstants, ChunkingConstants, BoundaryConfig


__all__ = ['calculate_single_shot_boundaries']


def calculate_single_shot_boundaries(
    bucket: str,
    key: str,
    filesize: int,
    total_chunks: int,
    config: BoundaryConfig,
    debug: bool = False
) -> Tuple[List[Tuple[int, int, bool]], int]:
    """
    Issue a single S3 GET at the middle of the file, derive a fixed window size,
    and return arithmetic shard boundaries.

    Args:
        bucket: S3 bucket name
        key: S3 object key
        filesize: Total file size in bytes
        total_chunks: Total number of chunks across all Lambdas
        config: BoundaryConfig (reads single_shot_* fields)
        debug: Enable debug output

    Returns:
        Tuple of (shard_ranges, window_size) where:
        - shard_ranges: List of (start_byte, end_byte, skip_first_line=False) tuples
        - window_size: Fixed window size in bytes (int) for Lambda correction
    """
    sample_size = config.single_shot_sample_kb * 1024
    safety_factor = config.single_shot_safety_factor

    if debug:
        print(f"[Single Shot] filesize={filesize}, total_chunks={total_chunks}, "
              f"sample_kb={config.single_shot_sample_kb}, safety_factor={safety_factor}")

    # Single sample at the middle of the file for representativeness
    avg_line_length = _sample_middle(
        bucket=bucket, key=key, filesize=filesize,
        sample_size=sample_size, debug=debug
    )

    # Derive window size, clamp to [MIN, MAX]
    window_size = max(
        int(avg_line_length * safety_factor),
        S3BoundaryConstants.MIN_WINDOW_SIZE_BYTES
    )
    window_size = min(window_size, S3BoundaryConstants.MAX_WINDOW_SIZE_BYTES)

    if debug:
        print(f"[Single Shot] avg_line_length={avg_line_length:.1f}B, "
              f"window_size={window_size}B ({window_size / ChunkingConstants.BYTES_TO_KB:.1f}KB)")

    # Arithmetic shard boundaries — identical pattern to adaptive_simple
    chunk_size = filesize // total_chunks
    shard_ranges = []

    for i in range(total_chunks):
        start = i * chunk_size
        end = filesize - 1 if i == total_chunks - 1 else (i + 1) * chunk_size - 1
        shard_ranges.append((start, end, False))  # Lambda will handle line alignment

    if debug:
        print(f"[Single Shot] Calculated {total_chunks} approximate boundaries, "
              f"Lambda will correct using {window_size / ChunkingConstants.BYTES_TO_KB:.1f}KB windows")

    return shard_ranges, window_size


def _sample_middle(bucket: str, key: str, filesize: int, sample_size: int, debug: bool) -> float:
    """
    Issue exactly one S3 GET at filesize // 2 and return the estimated average
    line length from that sample.

    Fallback logic (mirrors estimate_avg_line_length):
        - No newlines found  → 1 MB (FALLBACK_LONG_LINE_WINDOW_BYTES)
        - Any exception      → 64 KB (FALLBACK_WINDOW_SIZE_BYTES)
    """
    try:
        s3 = boto3.client('s3', region_name='us-east-1')

        start_pos = filesize // 2
        end_pos = min(start_pos + sample_size - 1, filesize - 1)

        if debug:
            print(f"[Sampling] File size: {filesize} bytes, taking 1 samples of {sample_size} bytes each")
            print(f"[Sampling] Sample 1 @ {start_pos} (middle of file)")

        response = s3.get_object(
            Bucket=bucket,
            Key=key,
            Range=f'bytes={start_pos}-{end_pos}'
        )
        sample = response['Body'].read()

        newline_count = sample.count(b'\n')

        if debug:
            sample_avg = len(sample) / newline_count if newline_count > 0 else 0
            print(f"[Sampling] Sample 1 @ {start_pos}: {newline_count} newlines, "
                  f"avg={sample_avg:.1f}B/line")

        if newline_count == 0:
            if debug:
                print(f"[Sampling] No newlines in {len(sample)} bytes across 1 samples")
                print(f"[Sampling] Defaulting to 1MB window (likely binary or very long lines)")
            return S3BoundaryConstants.FALLBACK_LONG_LINE_WINDOW_BYTES

        avg_line_length = len(sample) / newline_count

        if debug:
            print(f"[Sampling] Total: {len(sample)} bytes, {newline_count} newlines")
            print(f"[Sampling] Estimated avg line length: {avg_line_length:.1f} bytes")

        return avg_line_length

    except Exception as e:
        if debug:
            print(f"[Sampling] Error during sampling: {e}")
            print(f"[Sampling] Defaulting to 64KB window")
        return S3BoundaryConstants.FALLBACK_WINDOW_SIZE_BYTES
