"""
S3 Boundary Detection - Adaptive Simple Mode

Samples N equidistant positions from the file, computes mean average line length,
and derives a single fixed window size from that.  The window is handed to the
existing lambda-correction script (s3-chunk-reader-approx-correction.py)
as a plain integer — the same code path used by approx_with_correction.

Total S3 GET requests for boundary setup: N (default 5), regardless of the number
of shards.  Compare with the full adaptive mode which issues 1–2 GETs per shard.

Environment Variables:
    - USE_ADAPTIVE_SIMPLE=true: Enable this mode
    - PASH_ADAPTIVE_SIMPLE_NUM_SAMPLES: Number of equidistant sample points (default 5)
    - PASH_ADAPTIVE_SIMPLE_SAMPLE_KB: Bytes per sample in KB (default 256)
    - PASH_ADAPTIVE_SIMPLE_SAFETY_FACTOR: Multiplier on avg line length (default 1.5)

Functions:
    - calculate_adaptive_simple_boundaries: Sample → window → arithmetic shard_ranges

Usage:
    from compiler.serverless.s3_adaptive_simple import calculate_adaptive_simple_boundaries

    shard_ranges, window_size = calculate_adaptive_simple_boundaries(
        bucket="mybucket", key="file.txt", filesize=1000000,
        total_chunks=256, config=BoundaryConfig()
    )  # window_size is an int (bytes)
"""

import boto3
from typing import List, Tuple

from compiler.serverless.s3_config import S3BoundaryConstants, ChunkingConstants, BoundaryConfig
from compiler.serverless.s3_sampling import estimate_avg_line_length


__all__ = ['calculate_adaptive_simple_boundaries']


def calculate_adaptive_simple_boundaries(
    bucket: str,
    key: str,
    filesize: int,
    total_chunks: int,
    config: BoundaryConfig,
    debug: bool = False
) -> Tuple[List[Tuple[int, int, bool]], int]:
    """
    Sample the file at N equidistant points, derive a single fixed window size,
    and return arithmetic shard boundaries.

    Args:
        bucket: S3 bucket name
        key: S3 object key
        filesize: Total file size in bytes
        total_chunks: Total number of chunks across all Lambdas
        config: BoundaryConfig (reads adaptive_simple_* fields)
        debug: Enable debug output

    Returns:
        Tuple of (shard_ranges, window_size) where:
        - shard_ranges: List of (start_byte, end_byte, skip_first_line=False) tuples
        - window_size: Fixed window size in bytes (int) for Lambda correction
    """
    s3 = boto3.client('s3', region_name='us-east-1')

    num_samples = config.adaptive_simple_num_samples
    sample_size = config.adaptive_simple_sample_kb * 1024
    safety_factor = config.adaptive_simple_safety_factor

    if debug:
        print(f"[Adaptive Simple] filesize={filesize}, total_chunks={total_chunks}, "
              f"num_samples={num_samples}, sample_kb={config.adaptive_simple_sample_kb}, "
              f"safety_factor={safety_factor}")

    # Step 1: estimate average line length via shared sampling utility
    avg_line_length = estimate_avg_line_length(
        s3=s3, bucket=bucket, key=key, file_size=filesize,
        num_samples=num_samples, sample_size=sample_size,
        debug=debug
    )

    # Step 2: derive window size, clamp to [MIN, MAX]
    window_size = max(
        int(avg_line_length * safety_factor),
        S3BoundaryConstants.MIN_WINDOW_SIZE_BYTES
    )
    window_size = min(window_size, S3BoundaryConstants.MAX_WINDOW_SIZE_BYTES)

    if debug:
        print(f"[Adaptive Simple] avg_line_length={avg_line_length:.1f}B, "
              f"window_size={window_size}B ({window_size / ChunkingConstants.BYTES_TO_KB:.1f}KB)")

    # Step 3: arithmetic shard boundaries (same pattern as adaptive_dynamic / approx_correction)
    chunk_size = filesize // total_chunks
    shard_ranges = []

    for i in range(total_chunks):
        start = i * chunk_size
        end = filesize - 1 if i == total_chunks - 1 else (i + 1) * chunk_size - 1
        shard_ranges.append((start, end, False))  # Lambda will handle line alignment

    if debug:
        print(f"[Adaptive Simple] Calculated {total_chunks} approximate boundaries, "
              f"Lambda will correct using {window_size / ChunkingConstants.BYTES_TO_KB:.1f}KB windows")

    return shard_ranges, window_size
