#!/bin/bash

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

[ -z "$AWS_BUCKET" ] && {
  echo "AWS_BUCKET not set"
  exit
}

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/covid-mts"
INPUTS_DIR="$BENCHMARK_DIR/inputs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="covid-mts"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"

INPUTS=(
 # in.csv
 in_tiny.csv
  # in_small.csv
  # in_200M.csv
  # in_500M.csv
  # in_1G.csv
)

for INPUT in "${INPUTS[@]}"; do
  INPUT_PATH="$INPUTS_DIR/$INPUT"
  S3_URI=$S3_BUCKET_PREFIX/$S3_INPUTS_DIR/$INPUT
  echo "Uploading $INPUT_PATH to $S3_URI"
  aws s3 cp "$INPUT_PATH" "$S3_URI"
done
