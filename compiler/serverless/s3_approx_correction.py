"""
S3 Boundary Detection - Approximate with Lambda Correction Mode

This module implements the "approximate with lambda correction" mode where:
1. EC2 samples the file once to estimate average line length
2. EC2 calculates window size based on estimated line length
3. EC2 generates approximate (arithmetic) boundaries
4. Lambda workers use the window size to find exact line boundaries at runtime

This mode balances EC2 overhead (single sampling) with Lambda-side correction.

Environment Variables:
    - USE_APPROX_LAMBDA_CORRECTION=true: Enable this mode
    - PASH_BOUNDARY_WINDOW_MULTIPLIER: Target lines per window (default 500)

Functions:
    - calculate_approx_correction_boundaries: Sample file, calculate window, return approx boundaries

Usage:
    from compiler.serverless.s3_approx_correction import calculate_approx_correction_boundaries

    shard_ranges, window_size = calculate_approx_correction_boundaries(
        bucket="mybucket", key="file.txt", filesize=1000000,
        total_chunks=256, config=BoundaryConfig()
    )
"""

import boto3
from typing import List, Tuple

from compiler.serverless.s3_config import S3BoundaryConstants, ChunkingConstants, BoundaryConfig
from compiler.serverless.s3_sampling import estimate_avg_line_length
from compiler.serverless.ir_helper import time_block


__all__ = ['calculate_approx_correction_boundaries']


def calculate_approx_correction_boundaries(
    bucket: str,
    key: str,
    filesize: int,
    total_chunks: int,
    config: BoundaryConfig,
    debug: bool = False
) -> Tuple[List[Tuple[int, int, bool]], int]:
    """
    Sample file to estimate window size, then calculate approximate boundaries.

    Lambda workers will use the window_size to find exact line boundaries at runtime.

    Args:
        bucket: S3 bucket name
        key: S3 object key
        filesize: Total file size in bytes
        total_chunks: Total number of chunks across all Lambdas
        config: BoundaryConfig with window multiplier settings
        debug: Enable debug output

    Returns:
        Tuple of (shard_ranges, window_size) where:
        - shard_ranges: List of (start_byte, end_byte, skip_first_line=False) tuples
        - window_size: Window size in bytes for Lambda to use for correction
    """
    s3 = boto3.client('s3', region_name='us-east-1')

    # Sample file to estimate average line length
    with time_block("estimate_avg_line_length call"):
        avg_line_size = estimate_avg_line_length(
            s3=s3, bucket=bucket, key=key, file_size=filesize,
            num_samples=S3BoundaryConstants.DEFAULT_NUM_SAMPLES,
            sample_size=S3BoundaryConstants.SAMPLE_SIZE_BYTES,
            debug=debug
        )

    # Calculate window size based on target lines per window
    target_lines = config.boundary_window_multiplier
    window_size = int(avg_line_size * target_lines)
    window_size = max(
        S3BoundaryConstants.MIN_WINDOW_SIZE_BYTES,
        min(window_size, S3BoundaryConstants.MAX_WINDOW_SIZE_BYTES)
    )

    if debug:
        print(f"[Approx+Correction] Estimated avg line: {avg_line_size:.1f}B, "
              f"window: {window_size/ChunkingConstants.BYTES_TO_KB:.1f}KB")

    # Calculate approximate boundaries (pure arithmetic)
    chunk_size = filesize // total_chunks
    shard_ranges = []

    for i in range(total_chunks):
        start = i * chunk_size
        end = filesize - 1 if i == total_chunks - 1 else (i + 1) * chunk_size - 1
        shard_ranges.append((start, end, False))  # Lambda will handle line alignment

    if debug:
        print(f"[Approx+Correction] Calculated {total_chunks} approximate boundaries, "
              f"Lambda will correct using {window_size/ChunkingConstants.BYTES_TO_KB:.1f}KB windows")

    return shard_ranges, window_size
