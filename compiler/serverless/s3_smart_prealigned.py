"""
S3 Boundary Detection - Smart Boundaries Mode

This module implements the "smart boundaries" mode where the EC2 instance
pre-scans the S3 file to find exact line-aligned byte boundaries before
invoking Lambda workers. This minimizes Lambda-side overhead at the cost
of upfront EC2 scanning time.

Environment Variables:
    - USE_SMART_BOUNDARIES=true: Enable this mode
    - PASH_BOUNDARY_WINDOW_KB: Override auto-detected window size

Functions:
    - find_line_boundaries_smart: Download small windows around boundaries to find newlines
    - calculate_smart_boundaries: Wrapper that returns (shard_ranges, window_size=None)

Usage:
    from compiler.serverless.s3_smart_prealigned import calculate_smart_boundaries

    shard_ranges, window_size = calculate_smart_boundaries(
        bucket="mybucket", key="file.txt",
        total_lambdas=16, chunks_per_lambda=16
    )
"""

import os
import boto3
from typing import List, Tuple, Optional

from compiler.serverless.s3_config import S3BoundaryConstants, ChunkingConstants
from compiler.serverless.s3_sampling import estimate_avg_line_length, expand_search_for_newline
from compiler.serverless.ir_helper import time_block


__all__ = ['find_line_boundaries_smart', 'calculate_smart_boundaries']


def find_line_boundaries_smart(bucket, key, num_shards, chunks_per_lambda=1, window_size=None, debug=False):
    """
    Find line-aligned byte boundaries by downloading ONLY small windows around boundaries.

    If window_size is not provided, automatically determines optimal size by sampling file
    at multiple points (beginning, middle, end) to handle variable line lengths.

    Args:
        bucket: S3 bucket name
        key: S3 key
        num_shards: Number of lambda workers
        chunks_per_lambda: How many chunks each lambda processes
        window_size: Size of window to download around each boundary (None = auto-detect)
        debug: Enable debug output

    Returns:
        Tuple of (ranges, all_success) where:
        - ranges: List of (start_byte, end_byte, skip_first_line) tuples for each chunk
          - Has num_shards × chunks_per_lambda entries
          - skip_first_line=False if boundary was successfully line-aligned
          - skip_first_line=True if we had to fall back to approximate boundary
        - all_success: True if ALL boundaries were successfully found

    Example for 20GB file, 16 shards, 2 chunks per lambda (all successful):
        Total chunks: 16 × 2 = 32
        Downloads: 31 boundaries × 1MB = 31MB
        Returns: ([(0, 655359999, False), (655360000, 1310719999, False), ...], True)
    """
    with time_block("Total boundary scan"):
        s3 = boto3.client('s3')

        total_chunks = num_shards * chunks_per_lambda

        # Step 1: Get file size (needed for sampling positions)
        with time_block("S3 HEAD request"):
            response = s3.head_object(Bucket=bucket, Key=key)
            file_size = response['ContentLength']

        # Step 2: Adaptive window sizing
        if window_size is None:
            # Check for env var override first
            if 'PASH_BOUNDARY_WINDOW_KB' in os.environ:
                window_size_kb = int(os.environ.get('PASH_BOUNDARY_WINDOW_KB'))
                window_size = window_size_kb * ChunkingConstants.BYTES_TO_KB
                if debug:
                    print(f"[Boundary Scan] Using window_size={window_size/1024:.1f} KB (from env var)")
            else:
                # Multi-point sampling to estimate line length
                with time_block("Multi-point file sampling"):
                    # Take 5 samples of 256KB each = 1.25MB total
                    # Positions: 0%, 25%, 50%, 75%, 100% through file
                    avg_line_length = estimate_avg_line_length(
                        s3, bucket, key, file_size,
                        num_samples=S3BoundaryConstants.DEFAULT_NUM_SAMPLES,
                        sample_size=S3BoundaryConstants.SAMPLE_SIZE_BYTES,
                        debug=debug
                    )

                # Set window to hold ~500 lines (with 4KB minimum, 1MB maximum)
                # This provides safety margin while minimizing downloads
                target_lines_per_window = S3BoundaryConstants.DEFAULT_TARGET_LINES_PER_WINDOW
                window_size = int(avg_line_length * target_lines_per_window)
                window_size = max(S3BoundaryConstants.MIN_WINDOW_SIZE_BYTES,
                                min(window_size, S3BoundaryConstants.MAX_WINDOW_SIZE_BYTES))  # Clamp 4KB-1MB

                if debug:
                    print(f"[Boundary Scan] Adaptive window_size={window_size/1024:.1f} KB "
                          f"(avg_line={avg_line_length:.1f}B × {target_lines_per_window} lines)")

        all_boundaries_found = True  # Track if any boundary scan failed

        # Add verbosity level control
        verbose = debug and os.environ.get('PASH_VERBOSE_TIMING', 'false').lower() == 'true'

        if debug:
            print(f"[Boundary Scan] File size: {file_size} bytes, {num_shards} shards, {chunks_per_lambda} chunks/shard = {total_chunks} total chunks")

        # Edge case: Empty file
        if file_size == 0:
            if debug:
                print(f"[Boundary Scan] Empty file - returning trivial ranges")
            # All chunks get empty ranges, no skip needed
            return ([(0, 0, False) for _ in range(total_chunks)], True)

        # Edge case: Single chunk
        if total_chunks == 1:
            if debug:
                print(f"[Boundary Scan] Single chunk - no boundaries to scan")
            return ([(0, file_size - 1, False)], True)

        # Step 2: Calculate approximate chunk size
        approx_chunk_size = file_size // total_chunks

        # Track which boundaries were successfully found (vs fell back to approximate)
        # Index 0 is always 0 (start of file), so no skip needed for chunk 0
        boundary_success = [True]  # Chunk 0 starts at 0, always line-aligned

        boundaries = [0]  # First chunk starts at byte 0

        # Step 3: Find exact line boundaries
        with time_block(f"Boundary scan ({total_chunks-1} boundaries)"):
            for i in range(1, total_chunks):
                approx_boundary = i * approx_chunk_size
                found_newline = False

                # Download small window around the approximate boundary
                # 1MB window: 512KB before + 512KB after
                start = max(0, approx_boundary - window_size // 2)
                end = min(file_size - 1, approx_boundary + window_size // 2)

                if verbose:  # Only print per-chunk details if explicitly enabled
                    print(f"[Boundary Scan] Chunk {i}: downloading bytes {start}-{end} (window around {approx_boundary})")
                elif debug and i % 50 == 0:  # Print every 50th chunk as progress indicator
                    print(f"[Boundary Scan] Progress: {i}/{total_chunks} boundaries scanned")

                try:
                    response = s3.get_object(
                        Bucket=bucket,
                        Key=key,
                        Range=f'bytes={start}-{end}'
                    )
                    chunk = response['Body'].read()

                    # Find first newline AFTER the approximate boundary
                    offset_in_chunk = approx_boundary - start
                    try:
                        newline_pos = chunk.index(b'\n', offset_in_chunk)
                        # Boundary starts AFTER the newline (+1)
                        actual_boundary = start + newline_pos + 1
                        found_newline = True

                        if verbose:
                            print(f"[Boundary Scan] Chunk {i}: found newline at byte {actual_boundary}")

                        boundaries.append(actual_boundary)
                        boundary_success.append(True)

                    except ValueError:
                        # Newline not found in window - expand search
                        if verbose:
                            print(f"[Boundary Scan] Chunk {i}: newline not in window, expanding search...")

                        actual_boundary, found = expand_search_for_newline(
                            s3, bucket, key, approx_boundary, file_size
                        )
                        boundaries.append(actual_boundary)

                        if found:
                            boundary_success.append(True)
                            if verbose:
                                print(f"[Boundary Scan] Chunk {i}: found newline at byte {actual_boundary} (expanded search)")
                        else:
                            # No newline found even after expanding - use approximate, need skip
                            boundary_success.append(False)
                            all_boundaries_found = False
                            if verbose:
                                print(f"[Boundary Scan] Chunk {i}: NO NEWLINE FOUND, using approximate boundary {actual_boundary}, will skip first line")

                except Exception as e:
                    print(f"[Boundary Scan] ERROR at chunk {i}: {e}")
                    # Fallback to approximate boundary - will need skip_first_line
                    boundaries.append(approx_boundary)
                    boundary_success.append(False)
                    all_boundaries_found = False

        boundaries.append(file_size)  # Last chunk ends at EOF

        # Convert to (start, end, skip_first_line) ranges
        with time_block("Range calculation"):
            ranges = []
            for i in range(total_chunks):
                start_byte = boundaries[i]
                end_byte = boundaries[i + 1] - 1  # end_byte is inclusive

                # Validate range (handle edge case where boundaries collapse)
                if end_byte < start_byte:
                    # Empty or invalid range - give it at least one byte
                    end_byte = start_byte

                # skip_first_line is based on whether THIS chunk's START boundary was line-aligned
                # Chunk 0 starts at byte 0 (always line-aligned, no skip)
                # Chunk i>0 needs skip if boundary_success[i] is False
                skip_first_line = not boundary_success[i] if i > 0 else False

                ranges.append((start_byte, end_byte, skip_first_line))

        if debug:
            print(f"\n[Boundary Scan] Completed: {total_chunks} chunks, all_success={all_boundaries_found}")
            print(f"[Boundary Scan] Total data downloaded: ~{(total_chunks-1) * window_size / (1024*1024):.1f} MB")
            if verbose:  # Only print full range list if verbose mode enabled
                print(f"[Boundary Scan] Final ranges:")
                for i, (start, end, skip) in enumerate(ranges):
                    size_mb = (end - start + 1) / ChunkingConstants.BYTES_TO_MB
                    print(f"  Chunk {i}: bytes={start}-{end} ({size_mb:.2f} MB), skip={skip}")

        return (ranges, all_boundaries_found)


def calculate_smart_boundaries(bucket: str, key: str, total_lambdas: int,
                               chunks_per_lambda: int, debug: bool = False) -> Tuple[List[Tuple[int, int, bool]], None]:
    """
    Wrapper for smart boundaries mode - pre-scan file on EC2 to find exact line boundaries.

    Args:
        bucket: S3 bucket name
        key: S3 object key
        total_lambdas: Number of Lambda workers
        chunks_per_lambda: Chunks per Lambda
        debug: Enable debug output

    Returns:
        Tuple of (shard_ranges, window_size) where:
        - shard_ranges: List of (start_byte, end_byte, skip_first_line) tuples
        - window_size: None (not used by Lambda in smart mode)
    """
    shard_ranges, all_success = find_line_boundaries_smart(
        bucket=bucket,
        key=key,
        num_shards=total_lambdas,
        chunks_per_lambda=chunks_per_lambda,
        window_size=None,  # Auto-detect
        debug=debug
    )

    # Smart mode doesn't pass window_size to Lambda (boundaries are pre-aligned)
    return (shard_ranges, None)
