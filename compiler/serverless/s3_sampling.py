"""
S3 Boundary Detection - File Sampling Utilities

This module provides utilities for sampling S3 files to estimate characteristics
like average line length. These functions are shared by multiple boundary modes
(smart, approximate+correction).

Functions:
    - expand_search_for_newline: Exponentially expand search window to find newlines
    - estimate_avg_line_length: Sample file to estimate average line length

Usage:
    from compiler.serverless.s3_sampling import estimate_avg_line_length

    avg_len = estimate_avg_line_length(s3, bucket, key, file_size)
"""

import boto3, sys, os
from typing import Tuple

from serverless.s3_config import S3BoundaryConstants
from serverless.ir_helper import time_block


__all__ = ['expand_search_for_newline', 'estimate_avg_line_length']


def expand_search_for_newline(s3, bucket, key, start_pos, file_size, max_search=100*1024*1024):
    """
    Exponentially expand search window if newline not found in initial window.
    Handles edge cases like very long lines (>1MB).

    Args:
        s3: boto3 S3 client
        bucket: S3 bucket name
        key: S3 key
        start_pos: Starting byte position to search from
        file_size: Total file size
        max_search: Maximum bytes to search (default 100MB)

    Returns:
        Tuple of (position, found) where:
        - position: Byte position of newline +1, or end of search range if not found
        - found: True if newline was found, False otherwise
    """
    current_pos = start_pos
    chunk_size = S3BoundaryConstants.INITIAL_CHUNK_SIZE_BYTES  # Start with 1MB

    while current_pos < min(start_pos + max_search, file_size):
        try:
            end_pos = min(current_pos + chunk_size - 1, file_size - 1)
            response = s3.get_object(
                Bucket=bucket,
                Key=key,
                Range=f'bytes={current_pos}-{end_pos}'
            )
            chunk = response['Body'].read()

            newline_pos = chunk.index(b'\n')
            return (current_pos + newline_pos + 1, True)  # Found newline
        except ValueError:
            # No newline in this chunk, continue
            current_pos += chunk_size
            chunk_size = min(chunk_size * 2, S3BoundaryConstants.MAX_CHUNK_SIZE_BYTES)  # Double size, cap at 10MB
        except Exception as e:
            if 'InvalidRange' in str(e):
                break
            raise

    # No newline found in max_search range
    return (min(start_pos + max_search, file_size), False)  # Not found


def estimate_avg_line_length(s3, bucket, key, file_size, num_samples=S3BoundaryConstants.DEFAULT_NUM_SAMPLES,
                             sample_size=S3BoundaryConstants.SAMPLE_SIZE_BYTES, debug=False):
    """
    Estimate average line length by sampling multiple points across the file.

    Takes small samples from beginning, middle, end, and intermediate positions
    to get a representative average that handles variable line lengths.

    Args:
        s3: boto3 S3 client
        bucket: S3 bucket name
        key: S3 key
        file_size: Total file size in bytes
        num_samples: Number of sample points (default 5: start, 25%, 50%, 75%, end)
        sample_size: Bytes per sample (default 256KB)
        debug: Enable debug output

    Returns:
        Estimated average line length in bytes
    """
    try:
        total_bytes = 0
        total_newlines = 0

        # Calculate sample positions evenly distributed across file
        sample_positions = []
        if num_samples == 1:
            sample_positions = [0]
        else:
            for i in range(num_samples):
                # Distribute samples evenly: 0%, 25%, 50%, 75%, 100%
                position = int((file_size - sample_size) * i / (num_samples - 1))
                position = max(0, min(position, file_size - sample_size))
                sample_positions.append(position)

        if debug:
            print(f"[Sampling] File size: {file_size} bytes, taking {num_samples} samples of {sample_size} bytes each")

        for idx, start_pos in enumerate(sample_positions):
            end_pos = min(start_pos + sample_size - 1, file_size - 1)

            # Download sample
            response = s3.get_object(
                Bucket=bucket,
                Key=key,
                Range=f'bytes={start_pos}-{end_pos}'
            )
            sample = response['Body'].read()

            # Count newlines in this sample
            newline_count = sample.count(b'\n')
            total_newlines += newline_count
            total_bytes += len(sample)

            if debug:
                sample_avg = len(sample) / newline_count if newline_count > 0 else 0
                print(f"[Sampling] Sample {idx+1} @ {start_pos}: {newline_count} newlines, "
                      f"avg={sample_avg:.1f}B/line")

        if total_newlines == 0:
            # No newlines found - assume very long lines or binary data
            if debug:
                print(f"[Sampling] No newlines in {total_bytes} bytes across {num_samples} samples")
                print(f"[Sampling] Defaulting to 1MB window (likely binary or very long lines)")
            return S3BoundaryConstants.FALLBACK_LONG_LINE_WINDOW_BYTES  # Use 1MB window as safety

        avg_line_length = total_bytes / total_newlines

        if debug:
            print(f"[Sampling] Total: {total_bytes} bytes, {total_newlines} newlines")
            print(f"[Sampling] Estimated avg line length: {avg_line_length:.1f} bytes")

        return avg_line_length

    except Exception as e:
        if debug:
            print(f"[Sampling] Error during sampling: {e}")
            print(f"[Sampling] Defaulting to 64KB window")
        # Fallback to conservative default
        return S3BoundaryConstants.FALLBACK_WINDOW_SIZE_BYTES
