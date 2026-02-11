"""
S3 Gap Window Estimation - EC2-side utilities

This module estimates per-shard window sizes for the gap to the next newline
(residual life / inspection paradox aware) using byte-position sampling.
"""

import bisect
import random
from typing import List, Tuple
import sys, os

sys.path.append(os.path.join(os.getenv("PASH_TOP"), "compiler"))
from serverless.s3_config import S3BoundaryConstants


__all__ = [
    'sample_bytes',
    'estimate_gap_window',
    'compute_windows_per_shard',
]


def sample_bytes(s3_client, bucket: str, key: str, start: int, end: int) -> bytes:
    """Read a byte range from S3."""
    if end < start:
        return b""
    resp = s3_client.get_object(
        Bucket=bucket,
        Key=key,
        Range=f"bytes={int(start)}-{int(end)}"
    )
    return resp['Body'].read()


def _find_newline_positions(buf: bytes) -> List[int]:
    return [i for i, b in enumerate(buf) if b == 10]


def _max_run_length(buf: bytes) -> int:
    if not buf:
        return 0
    return max((len(seg) for seg in buf.split(b"\n")), default=0)


def estimate_gap_window(
    buf: bytes,
    delta: float,
    k_samples: int,
    safety_factor: float,
    clamp_min: int,
    clamp_max: int,
    debug: bool = False,
) -> int:
    """
    Estimate a window size for the gap to the next newline.

    Args:
        buf: Sample buffer.
        delta: Target tail probability (P(G > w) <= delta).
        k_samples: Number of random byte offsets to sample.
        safety_factor: Multiplier for max_run safety floor.
        clamp_min: Minimum window size in bytes.
        clamp_max: Maximum window size in bytes.
    """
    if not buf:
        return clamp_min

    newline_positions = _find_newline_positions(buf)
    if not newline_positions:
        if debug:
            print("[Gap Window] No newlines in sample; using clamp_max")
        return clamp_max

    last_nl = newline_positions[-1]
    if last_nl <= 0:
        return max(clamp_min, 1)

    # Clamp delta to [0, 1)
    if delta < 0.0:
        delta = 0.0
    if delta >= 1.0:
        delta = 0.999999

    gaps = []
    for _ in range(k_samples):
        x = random.randint(0, last_nl)
        idx = bisect.bisect_left(newline_positions, x)
        if idx >= len(newline_positions):
            idx = len(newline_positions) - 1
        gap = newline_positions[idx] - x + 1
        gaps.append(gap)

    gaps.sort()
    q_index = int((1.0 - delta) * (len(gaps) - 1))
    q_index = max(0, min(q_index, len(gaps) - 1))
    q = gaps[q_index]

    max_run = _max_run_length(buf)
    safety = int(max_run * safety_factor)

    window = max(q, safety, clamp_min)
    window = min(window, clamp_max)

    if debug:
        print(
            f"[Gap Window] q={q} safety={safety} "
            f"clamp=[{clamp_min},{clamp_max}] -> window={window}"
        )

    return window


def _choose_sample_centers(start: int, end: int, sample_size: int, max_centers: int = 2) -> List[int]:
    shard_len = max(0, end - start + 1)
    if shard_len <= sample_size or max_centers <= 1:
        return [start + shard_len // 2]

    if shard_len <= sample_size * 2:
        return [start + shard_len // 2]

    centers = [start + shard_len // 4, start + (3 * shard_len) // 4]
    return centers[:max_centers]


def compute_windows_per_shard(
    s3_client,
    bucket: str,
    key: str,
    shard_ranges: List[Tuple[int, int]],
    sample_kb: int,
    delta: float,
    k_samples: int,
    safety_factor: float,
    max_window_kb: int,
    debug: bool = False,
) -> List[int]:
    """
    Compute a per-shard window vector for gap-to-next-newline.

    Each shard gets 1-2 samples; the shard window is the max of sample windows.
    """
    sample_size = int(sample_kb) * 1024
    clamp_min = S3BoundaryConstants.MIN_WINDOW_SIZE_BYTES
    clamp_max = int(max_window_kb) * 1024

    windows = []
    for shard_index, (start, end) in enumerate(shard_ranges):
        centers = _choose_sample_centers(start, end, sample_size, max_centers=2)
        shard_windows = []

        for center in centers:
            sample_start = max(start, center - sample_size // 2)
            sample_end = min(end, sample_start + sample_size - 1)
            buf = sample_bytes(s3_client, bucket, key, sample_start, sample_end)
            w = estimate_gap_window(
                buf,
                delta=delta,
                k_samples=k_samples,
                safety_factor=safety_factor,
                clamp_min=clamp_min,
                clamp_max=clamp_max,
                debug=debug,
            )
            shard_windows.append(w)

        window = max(shard_windows) if shard_windows else clamp_max
        windows.append(window)

        if debug:
            print(f"[Gap Window] Shard {shard_index}: window={window} bytes")

    return windows
