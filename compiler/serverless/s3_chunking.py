"""
S3 Boundary Detection - Chunk Distribution and Formatting

This module handles distribution of byte range chunks to Lambda workers and
formatting of byte range parameters. Used by all boundary modes.

Functions:
    - get_s3_size: Get S3 object size
    - distribute_chunks_to_lambda: Round-robin chunk distribution to lambda workers
    - format_byte_range_parameter: Format byte ranges for Lambda invocation

Usage:
    from compiler.serverless.s3_chunking import distribute_chunks_to_lambda

    chunks = distribute_chunks_to_lambda(0, 16, 16, shard_ranges)
"""

import json
import boto3
from typing import Dict, List, Tuple, Any

from serverless.s3_config import ChunkingConstants


__all__ = ['get_s3_size', 'distribute_chunks_to_lambda', 'format_byte_range_parameter']


def get_s3_size(bucket, key):
    """Get S3 object size."""
    s3 = boto3.client('s3')
    response = s3.head_object(Bucket=bucket, Key=key)
    return response['ContentLength']


def distribute_chunks_to_lambda(
    lambda_counter: int,
    total_lambdas: int,
    chunks_per_lambda: int,
    shard_ranges: List[Tuple[int, int, bool]],
    debug: bool = False
) -> List[Dict[str, Any]]:
    """
    Distribute chunks to a single lambda using round-robin (interleaved) strategy.

    Formula: global_chunk_id = lambda_counter + (chunk_index * total_lambdas)

    Example (16 lambdas, 256 total chunks, 16 chunks/lambda):
        Lambda 0: blocks [0, 16, 32, 48, ..., 240]
        Lambda 1: blocks [1, 17, 33, 49, ..., 241]

    Args:
        lambda_counter: Current lambda index (0-based)
        total_lambdas: Total number of lambdas
        chunks_per_lambda: How many chunks each lambda processes
        shard_ranges: List of (start_byte, end_byte, skip_first_line) tuples
        debug: Enable debug logging

    Returns:
        List of chunk dicts with keys: start, end, block_id, [skip_first_line]
    """
    lambda_chunks = []

    for chunk_index in range(chunks_per_lambda):
        # CRITICAL: Round-robin distribution formula
        global_chunk_id = lambda_counter + (chunk_index * total_lambdas)

        # Safety: handle edge case where chunks don't divide evenly
        if global_chunk_id >= len(shard_ranges):
            if debug:
                print(f"[Chunk Dist] Lambda {lambda_counter} chunk {chunk_index}: "
                      f"ID {global_chunk_id} >= {len(shard_ranges)}, skipping")
            break

        start_byte, end_byte, skip_first_line = shard_ranges[global_chunk_id]

        chunk_dict = {
            "start": start_byte,
            "end": end_byte,
            "block_id": global_chunk_id,  # CRITICAL: Sequential ID for r_merge
            "shard_id": chunk_index  # Coarse shard index for per-shard window stats
        }

        # Only include skip_first_line if True (smart boundaries mode)
        if skip_first_line:
            chunk_dict["skip_first_line"] = skip_first_line

        lambda_chunks.append(chunk_dict)

    if debug and lambda_chunks:
        chunk_ids = [c['block_id'] for c in lambda_chunks]
        total_mb = sum((c['end'] - c['start'] + 1) / ChunkingConstants.BYTES_TO_MB
                      for c in lambda_chunks)
        print(f"[Chunk Dist] Lambda {lambda_counter}: {len(lambda_chunks)} chunks, "
              f"IDs {chunk_ids}, {total_mb:.2f} MB total")

    return lambda_chunks


def format_byte_range_parameter(
    lambda_chunks: List[Dict[str, Any]],
    chunks_per_lambda: int,
    lambda_counter: int = 0,
    debug: bool = False,
    force_json: bool = False,
) -> str:
    """
    Format byte range parameter for Lambda invocation.

    Single-chunk mode: Returns "bytes=START-END" (unless force_json=True)
    Multi-chunk mode: Returns JSON array string (single-quoted for shell safety)

    Args:
        lambda_chunks: List of chunk dicts (start, end, block_id, skip_first_line)
        chunks_per_lambda: Number of chunks per lambda (determines format)
        lambda_counter: Lambda index (for debug logging)
        debug: Enable debug logging
        force_json: Force JSON output even for single-chunk mode

    Returns:
        Formatted byte_range string ready for Lambda invocation
    """
    if chunks_per_lambda == 1 and not force_json:
        # Single-chunk mode: backward compatible string format
        chunk = lambda_chunks[0]
        byte_range = f"bytes={chunk['start']}-{chunk['end']}"

        if debug:
            skip_info = f", skip={chunk.get('skip_first_line', False)}" \
                if 'skip_first_line' in chunk else ""
            print(f"[Byte Range] Lambda {lambda_counter}: {byte_range}{skip_info}")

        return byte_range
    else:
        # Multi-chunk (or forced JSON) mode: JSON array (shell-safe with single quotes)
        byte_range_json = json.dumps(lambda_chunks)
        byte_range = f"'{byte_range_json}'"

        if debug:
            total_mb = sum((c['end'] - c['start'] + 1) / ChunkingConstants.BYTES_TO_MB
                          for c in lambda_chunks)
            chunk_ids = [c['block_id'] for c in lambda_chunks]
            print(f"[Byte Range] Lambda {lambda_counter}: {len(lambda_chunks)} chunks, "
                  f"IDs {chunk_ids}, {total_mb:.2f} MB total")

        return byte_range
