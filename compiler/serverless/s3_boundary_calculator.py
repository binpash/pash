"""
S3 Boundary Detection - Facade Module

This module provides a unified interface (BoundaryCalculator) that dispatches
to the appropriate boundary mode implementation based on configuration.

The BoundaryCalculator class:
- Reads configuration to determine which boundary mode to use
- Dispatches to mode-specific modules
- Caches results across multiple Lambda invocations
- Provides a clean API for ir_helper.py

Boundary Modes (priority order):
1. Adaptive: Quantile analysis with sampling (USE_ADAPTIVE_BOUNDARIES=true)
2. Dynamic: Exponential expansion (USE_DYNAMIC_BOUNDARIES=true)
3. Adaptive Simple: N samples → single fixed window (USE_ADAPTIVE_SIMPLE=true)
4. Single Shot: 1 sample at midpoint, higher safety factor (USE_SINGLE_SHOT=true)
5. Approx+Correction: Sample once, Lambda corrects (USE_APPROX_LAMBDA_CORRECTION=true)
6. Smart: EC2 pre-scans boundaries (USE_SMART_BOUNDARIES=true)
7. Legacy Approximate: No pre-calculation (default)

Usage:
    from compiler.serverless.s3_boundary_calculator import BoundaryCalculator
    from compiler.serverless.s3_config import BoundaryConfig

    config = BoundaryConfig()
    calculator = BoundaryCalculator(config, debug=True)
    shard_ranges, window_size = calculator.calculate_boundaries(
        bucket="mybucket", key="file.txt", filesize=1000000,
        total_lambdas=16, chunks_per_lambda=16
    )
"""

from typing import List, Tuple, Any, Optional

from serverless.s3_config import BoundaryConfig


__all__ = ['BoundaryCalculator']


class BoundaryCalculator:
    """Calculate S3 byte boundaries once per job (cached across lambdas)."""

    def __init__(self, config: BoundaryConfig, debug: bool = False):
        self.config = config
        self.debug = debug
        self.shard_ranges = None
        self.cached_window_size = None

    def calculate_boundaries(
        self,
        bucket: str,
        key: str,
        filesize: int,
        total_lambdas: int,
        chunks_per_lambda: int
    ) -> Tuple[List[Tuple[int, int, bool]], Any]:
        """
        Calculate boundaries once per job.

        Dispatches to the appropriate boundary mode based on config.
        Results are cached so subsequent calls return immediately.

        Args:
            bucket: S3 bucket name
            key: S3 object key
            filesize: Total file size in bytes
            total_lambdas: Number of Lambda workers
            chunks_per_lambda: Chunks per Lambda worker

        Returns:
            Tuple of (shard_ranges, window_size_param) where:
            - shard_ranges: List of (start_byte, end_byte, skip_first_line) tuples
            - window_size_param: Mode-specific parameter for Lambda
                - vector string for adaptive mode
                - None for dynamic/smart mode
                - int (bytes) for approx+correction mode
        """
        if self.shard_ranges is not None:
            # Already computed - reuse cached result
            return self.shard_ranges, self.cached_window_size

        total_chunks = total_lambdas * chunks_per_lambda

        # Dispatch to appropriate mode (priority order matches BoundaryConfig)
        if self.config.use_adaptive_boundaries or self.config.use_dynamic_boundaries:
            # Adaptive or Dynamic mode: Pure arithmetic boundaries
            from compiler.serverless.s3_adaptive_dynamic import calculate_adaptive_dynamic_boundaries

            self.shard_ranges, self.cached_window_size = calculate_adaptive_dynamic_boundaries(
                use_adaptive=self.config.use_adaptive_boundaries,
                filesize=filesize,
                total_chunks=total_chunks,
                bucket=bucket,
                key=key,
                total_lambdas=total_lambdas,
                chunks_per_lambda=chunks_per_lambda,
                config=self.config,
                debug=self.debug
            )

        elif self.config.use_adaptive_simple:
            # Adaptive-simple mode: N samples → single fixed window → lambda correction
            from compiler.serverless.s3_adaptive_simple import calculate_adaptive_simple_boundaries

            self.shard_ranges, self.cached_window_size = calculate_adaptive_simple_boundaries(
                bucket=bucket,
                key=key,
                filesize=filesize,
                total_chunks=total_chunks,
                config=self.config,
                debug=self.debug
            )

        elif self.config.use_single_shot:
            # Single-shot mode: 1 sample at file midpoint → fixed window → lambda correction
            from compiler.serverless.s3_single_shot import calculate_single_shot_boundaries

            self.shard_ranges, self.cached_window_size = calculate_single_shot_boundaries(
                bucket=bucket,
                key=key,
                filesize=filesize,
                total_chunks=total_chunks,
                config=self.config,
                debug=self.debug
            )

        elif self.config.use_approx_with_correction:
            # Approximate + Correction mode: Sample file, calculate window, return approx boundaries
            from compiler.serverless.s3_approx_correction import calculate_approx_correction_boundaries

            self.shard_ranges, self.cached_window_size = calculate_approx_correction_boundaries(
                bucket=bucket,
                key=key,
                filesize=filesize,
                total_chunks=total_chunks,
                config=self.config,
                debug=self.debug
            )

        elif self.config.use_smart_boundaries:
            # Smart mode: EC2-side boundary scanning
            from compiler.serverless.s3_smart_prealigned import calculate_smart_boundaries

            self.shard_ranges, self.cached_window_size = calculate_smart_boundaries(
                bucket=bucket,
                key=key,
                total_lambdas=total_lambdas,
                chunks_per_lambda=chunks_per_lambda,
                debug=self.debug
            )

        else:
            # Legacy approximate mode - no pre-calculation
            # Boundaries will be calculated inline in ir_helper.py
            self.shard_ranges = None
            self.cached_window_size = None

        return self.shard_ranges, self.cached_window_size
