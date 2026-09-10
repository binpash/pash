#!/bin/bash
# Minimal single serverless run for the agent benchmark (bypasses the matrix).
#   PASH_TOP and AWS_BUCKET must be set; inputs must already be on S3 (./upload.sh).
cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ]  && { echo "PASH_TOP not set";  exit 1; }
[ -z "$AWS_BUCKET" ] && { echo "AWS_BUCKET not set"; exit 1; }

WIDTH=${WIDTH:-1}
IN="agent/inputs/logs/"                 # S3 key prefix
OUT="agent/outputs/example.sh:logs:${WIDTH}:leash-"

echo "IN=$IN OUT=$OUT $PASH_TOP/pa.sh --serverless_exec -w$WIDTH scripts/example.sh"
time IN="$IN" OUT="$OUT" "$PASH_TOP/pa.sh" --serverless_exec -w"$WIDTH" scripts/example.sh

sleep 15
logs_dir="logs/example.sh:logs:${WIDTH}"
rm -rf "$logs_dir"
python3 "$PASH_TOP/scripts/serverless/utils.py" "$logs_dir" || true
