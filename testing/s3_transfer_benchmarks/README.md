# S3 Transfer Benchmarks

This folder contains the organized benchmark assets used to compare data transfer paths for serverless PaSh testing.

## What is measured

1. `S3 -> EC2` full object download timing
2. `S3 -> Lambda` full object download timing
3. `S3 -> Lambda` byte-range download timing
4. `S3 -> EC2 -> Lambda` streaming timing

## Layout

- `ec2_s3_benchmark.py`: direct S3-to-EC2 pull benchmark script
- `orchestrators/`: manual orchestrators used to invoke Lambda sort workers
- `lambda_workers/`: Lambda worker handlers and deployment scripts
- `analysis/plot.py`: plotting script for benchmark results
- `analysis/*.png`: generated comparison figures

## Quick usage

From this directory:

```bash
# Deploy worker lambdas
./lambda_workers/deploy_lambda_sort.sh
./lambda_workers/deploy_lambda_sort_byte_ranges.sh

# Run byte-range orchestrator
python3 orchestrators/manual_s3_orchestrator_byte_ranges.py \
  --bucket "$AWS_BUCKET" \
  --input oneliners/inputs/1G.txt \
  --output oneliners/outputs/byte-range-result.txt \
  --workers 2

# Plot current benchmark summary
python3 analysis/plot.py
```

## Notes

- Legacy/ad-hoc experimental artifacts were moved under `testing/legacy/` and are intentionally not committed.
- Deployment scripts package Lambda zips locally in `lambda_workers/`.
