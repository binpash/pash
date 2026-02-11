"""
S3 Boundary Detection - Adaptive and Dynamic Modes

This module implements both adaptive and dynamic boundary modes. These modes share
identical EC2-side logic (pure arithmetic boundaries) and only differ in what sentinel
value is passed to Lambda workers.

EC2-Side (this module):
    - Both modes: Calculate approximate boundaries using filesize / total_chunks
    - Adaptive mode: Computes shard-specific gap windows and returns vector string
    - Dynamic mode: Returns window_size=None

Lambda-Side (different behaviors):
    - Adaptive mode: Use shard-specific overlap windows (no lambda-side sampling)
    - Dynamic mode: Use exponential expansion starting from 1KB

Environment Variables:
    - USE_ADAPTIVE_BOUNDARIES=true: Enable adaptive mode
    - USE_DYNAMIC_BOUNDARIES=true: Enable dynamic mode

Functions:
    - calculate_adaptive_dynamic_boundaries: Shared EC2-side logic for both modes

Usage:
    from compiler.serverless.s3_adaptive_dynamic import calculate_adaptive_dynamic_boundaries

    # Adaptive mode
    shard_ranges, window_size = calculate_adaptive_dynamic_boundaries(
        use_adaptive=True, filesize=1000000, total_chunks=256,
        bucket="mybucket", key="file.txt", total_lambdas=16, chunks_per_lambda=16,
        config=BoundaryConfig()
    )  # window_size = "w0,w1,..."

    # Dynamic mode
    shard_ranges, window_size = calculate_adaptive_dynamic_boundaries(
        use_adaptive=False, filesize=1000000, total_chunks=256
    )  # window_size = None
"""

from typing import List, Tuple, Optional


__all__ = ['calculate_adaptive_dynamic_boundaries']


def calculate_adaptive_dynamic_boundaries(
    use_adaptive: bool,
    filesize: int,
    total_chunks: int,
    debug: bool = False,
    *,
    bucket: Optional[str] = None,
    key: Optional[str] = None,
    total_lambdas: Optional[int] = None,
    chunks_per_lambda: Optional[int] = None,
    config=None,
) -> Tuple[List[Tuple[int, int, bool]], Optional[str]]:
    """
    Calculate approximate boundaries for adaptive or dynamic mode.

    Both modes use identical EC2-side logic (pure arithmetic). The only difference
    is the sentinel value passed to Lambda workers.

    Args:
        use_adaptive: True for adaptive mode, False for dynamic mode
        filesize: Total file size in bytes
        total_chunks: Total number of chunks across all Lambdas
        debug: Enable debug output

    Returns:
        Tuple of (shard_ranges, window_size) where:
        - shard_ranges: List of (start_byte, end_byte, skip_first_line=False) tuples
        - window_size: vector string for adaptive mode, None for dynamic mode
    """
    # Pure arithmetic boundaries - no S3 access needed
    chunk_size = filesize // total_chunks
    shard_ranges = []

    for i in range(total_chunks):
        start = i * chunk_size
        end = filesize - 1 if i == total_chunks - 1 else (i + 1) * chunk_size - 1
        shard_ranges.append((start, end, False))  # Lambda will handle line alignment

    window_size = None
    mode_name = "adaptive" if use_adaptive else "dynamic"

    if use_adaptive:
        if bucket is None or key is None or chunks_per_lambda is None or config is None:
            # Fallback to legacy adaptive sentinel if required inputs are missing
            window_size = "adaptive"
            if debug:
                print("[Adaptive Mode] Missing inputs for gap windows; using legacy sentinel")
            return shard_ranges, window_size

        import boto3, sys, os
        sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
        from serverless.s3_gap_quantile_windows import compute_windows_per_shard

        # Build coarse shard ranges (one per chunk_index / shard_id)
        shard_count = max(1, int(chunks_per_lambda))
        shard_size = filesize // shard_count
        coarse_ranges = []
        for i in range(shard_count):
            start = i * shard_size
            end = filesize - 1 if i == shard_count - 1 else (i + 1) * shard_size - 1
            coarse_ranges.append((start, end))

        s3_client = boto3.client('s3', region_name='us-east-1')
        windows = compute_windows_per_shard(
            s3_client=s3_client,
            bucket=bucket,
            key=key,
            shard_ranges=coarse_ranges,
            sample_kb=config.gap_sample_kb,
            delta=config.gap_delta,
            k_samples=config.gap_k_samples,
            safety_factor=config.gap_safety_factor,
            max_window_kb=config.gap_max_window_kb,
            debug=debug,
        )

        window_size = ",".join(str(w) for w in windows)

    if debug:
        print(f"[{mode_name.capitalize()} Mode] Calculated {total_chunks} approximate boundaries")
        if use_adaptive:
            print(f"[{mode_name.capitalize()} Mode] Lambda will use shard-specific gap windows")
        else:
            print(f"[{mode_name.capitalize()} Mode] Lambda will use {mode_name} boundary detection")

    return shard_ranges, window_size
