"""
S3 Boundary Detection - Configuration and Constants

This module contains shared configuration and constants used across all
S3 boundary detection modes (smart, adaptive, dynamic, approximate).

Classes:
    - BoundaryConfig: Reads boundary mode configuration from environment variables
    - S3BoundaryConstants: Constants for S3 boundary scanning and file sampling
    - ChunkingConstants: Constants for chunk distribution

Usage:
    from compiler.serverless.s3_config import BoundaryConfig, S3BoundaryConstants

    config = BoundaryConfig()
    if config.use_smart_boundaries:
        # Use smart mode
        pass
"""

import os
from typing import Optional


__all__ = ['BoundaryConfig', 'S3BoundaryConstants', 'ChunkingConstants']


class S3BoundaryConstants:
    """Constants for S3 boundary scanning and file sampling."""
    # File sampling
    DEFAULT_NUM_SAMPLES = 4
    DEFAULT_SAMPLE_SIZE_KB = 256
    SAMPLE_SIZE_BYTES = DEFAULT_SAMPLE_SIZE_KB * 1024

    # Window sizing
    MIN_WINDOW_SIZE_KB = 4
    MAX_WINDOW_SIZE_KB = 1024
    MIN_WINDOW_SIZE_BYTES = MIN_WINDOW_SIZE_KB * 1024
    MAX_WINDOW_SIZE_BYTES = MAX_WINDOW_SIZE_KB * 1024
    DEFAULT_TARGET_LINES_PER_WINDOW = 500

    # Gap-window estimation (adaptive mode v2)
    DEFAULT_GAP_SAMPLE_KB = 256
    DEFAULT_GAP_DELTA = 1e-3
    DEFAULT_GAP_K_SAMPLES = 4096
    DEFAULT_GAP_SAFETY_FACTOR = 1.2
    DEFAULT_GAP_MAX_WINDOW_KB = 1024

    # Search limits
    MAX_NEWLINE_SEARCH_MB = 100
    MAX_NEWLINE_SEARCH_BYTES = MAX_NEWLINE_SEARCH_MB * 1024 * 1024
    INITIAL_CHUNK_SIZE_MB = 1
    INITIAL_CHUNK_SIZE_BYTES = INITIAL_CHUNK_SIZE_MB * 1024 * 1024
    MAX_CHUNK_SIZE_MB = 10
    MAX_CHUNK_SIZE_BYTES = MAX_CHUNK_SIZE_MB * 1024 * 1024

    # Fallback values
    FALLBACK_WINDOW_SIZE_KB = 64
    FALLBACK_WINDOW_SIZE_BYTES = FALLBACK_WINDOW_SIZE_KB * 1024
    FALLBACK_LONG_LINE_WINDOW_MB = 1
    FALLBACK_LONG_LINE_WINDOW_BYTES = FALLBACK_LONG_LINE_WINDOW_MB * 1024 * 1024


class ChunkingConstants:
    """Constants for chunk distribution."""
    DEFAULT_CHUNKS_PER_LAMBDA = 1
    BYTES_TO_MB = 1024 * 1024
    BYTES_TO_KB = 1024


class BoundaryConfig:
    """Read and cache boundary mode configuration from environment."""

    def __init__(self):
        # Debug flags
        self._debug_enabled = os.environ.get('PASH_DEBUG', 'false').lower() == 'true'
        self._verbose_enabled = self._debug_enabled and \
            os.environ.get('PASH_VERBOSE_TIMING', 'false').lower() == 'true'

        # Boundary mode selection (priority order)
        self.use_adaptive_boundaries = \
            os.environ.get('USE_ADAPTIVE_BOUNDARIES', 'false').lower() == 'true'
        self.use_dynamic_boundaries = \
            os.environ.get('USE_DYNAMIC_BOUNDARIES', 'false').lower() == 'true'
        self.use_adaptive_simple = \
            os.environ.get('USE_ADAPTIVE_SIMPLE', 'false').lower() == 'true'
        self.use_single_shot = \
            os.environ.get('USE_SINGLE_SHOT', 'false').lower() == 'true'
        self.use_approx_with_correction = \
            os.environ.get('USE_APPROX_LAMBDA_CORRECTION', 'false').lower() == 'true'
        self.use_smart_boundaries = \
            os.environ.get('USE_SMART_BOUNDARIES', 'false').lower() == 'true'

        # Chunk configuration
        self.chunks_per_lambda = self._read_int('PASH_S3_CHUNKS_PER_LAMBDA',
                                                ChunkingConstants.DEFAULT_CHUNKS_PER_LAMBDA,
                                                min_val=1)
        self.aws_bucket = os.environ.get('AWS_BUCKET')

        # Adaptive mode parameters
        self.adaptive_target_lines = self._read_int('PASH_ADAPTIVE_TARGET_LINES',
                                                    S3BoundaryConstants.DEFAULT_TARGET_LINES_PER_WINDOW)
        self.adaptive_retry_risk = self._read_float('PASH_ADAPTIVE_RETRY_RISK', 0.001)
        self.adaptive_sample_kb = self._read_int('PASH_ADAPTIVE_SAMPLE_KB',
                                                S3BoundaryConstants.DEFAULT_SAMPLE_SIZE_KB)

        # Adaptive-simple mode parameters
        self.adaptive_simple_num_samples = self._read_int('PASH_ADAPTIVE_SIMPLE_NUM_SAMPLES',
                                                           S3BoundaryConstants.DEFAULT_NUM_SAMPLES,
                                                           min_val=1)
        self.adaptive_simple_sample_kb = self._read_int('PASH_ADAPTIVE_SIMPLE_SAMPLE_KB',
                                                         S3BoundaryConstants.DEFAULT_SAMPLE_SIZE_KB,
                                                         min_val=1)
        self.adaptive_simple_safety_factor = self._read_float('PASH_ADAPTIVE_SIMPLE_SAFETY_FACTOR', 1.5)

        # Single-shot mode parameters (1 sample at file midpoint)
        self.single_shot_sample_kb = self._read_int('PASH_SINGLE_SHOT_SAMPLE_KB',
                                                     S3BoundaryConstants.DEFAULT_SAMPLE_SIZE_KB,
                                                     min_val=1)
        self.single_shot_safety_factor = self._read_float('PASH_SINGLE_SHOT_SAFETY_FACTOR', 2.0)

        # Gap-window adaptive parameters (EC2-side)
        self.gap_sample_kb = self._read_int('PASH_GAP_SAMPLE_KB',
                                            S3BoundaryConstants.DEFAULT_GAP_SAMPLE_KB,
                                            min_val=1)
        self.gap_delta = self._read_float('PASH_GAP_DELTA',
                                          S3BoundaryConstants.DEFAULT_GAP_DELTA)
        self.gap_k_samples = self._read_int('PASH_GAP_K_SAMPLES',
                                            S3BoundaryConstants.DEFAULT_GAP_K_SAMPLES,
                                            min_val=1)
        self.gap_safety_factor = self._read_float('PASH_GAP_SAFETY_FACTOR',
                                                  S3BoundaryConstants.DEFAULT_GAP_SAFETY_FACTOR)
        self.gap_max_window_kb = self._read_int('PASH_GAP_MAX_WINDOW_KB',
                                                S3BoundaryConstants.DEFAULT_GAP_MAX_WINDOW_KB,
                                                min_val=1)

        # Window size configuration
        self.boundary_window_kb = self._read_int_optional('PASH_BOUNDARY_WINDOW_KB')
        self.boundary_window_multiplier = self._read_int('PASH_BOUNDARY_WINDOW_MULTIPLIER',
                                                         S3BoundaryConstants.DEFAULT_TARGET_LINES_PER_WINDOW)

    @staticmethod
    def _read_int(key: str, default: int, min_val: int = None) -> int:
        """Read integer from env with validation."""
        value = int(os.environ.get(key, str(default)))
        if min_val is not None and value < min_val:
            print(f"[Config] WARNING: {key}={value} invalid, using {default}")
            return default
        return value

    @staticmethod
    def _read_int_optional(key: str) -> Optional[int]:
        """Read optional integer from env."""
        value = os.environ.get(key)
        return int(value) if value else None

    @staticmethod
    def _read_float(key: str, default: float) -> float:
        """Read float from env."""
        return float(os.environ.get(key, str(default)))

    def is_debug(self) -> bool:
        return self._debug_enabled

    def is_verbose(self) -> bool:
        return self._verbose_enabled

    def get_boundary_mode_name(self) -> str:
        """Get human-readable name of active boundary mode."""
        if self.use_adaptive_boundaries:
            return "adaptive"
        elif self.use_dynamic_boundaries:
            return "dynamic"
        elif self.use_adaptive_simple:
            return "adaptive (simple)"
        elif self.use_single_shot:
            return "single shot (1 sample)"
        elif self.use_approx_with_correction:
            return "approximate with lambda correction"
        elif self.use_smart_boundaries:
            return "smart (pre-aligned)"
        else:
            return "approximate (legacy)"
