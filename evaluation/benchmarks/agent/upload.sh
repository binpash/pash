#!/bin/bash
# Upload the log corpus to S3 so the serverless (Lambda) run can read it.
# Mirrors oneliners/upload.sh and nlp/upload.sh.
cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && { echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"; exit 1; }
[ -z "$AWS_BUCKET" ] && { echo "AWS_BUCKET not set"; exit 1; }

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/agent"
LOGS_DIR="$BENCHMARK_DIR/inputs/logs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="agent"
S3_LOGS_DIR="$S3_BENCHMARK_DIR/inputs/logs"

[ -d "$LOGS_DIR" ] || { echo "no inputs at $LOGS_DIR; run ./inputs.sh first"; exit 1; }

echo "Uploading $LOGS_DIR/ to $S3_BUCKET_PREFIX/$S3_LOGS_DIR/"
aws s3 sync "$LOGS_DIR/" "$S3_BUCKET_PREFIX/$S3_LOGS_DIR/" --exclude '*' --include '*.log'
echo "done"
