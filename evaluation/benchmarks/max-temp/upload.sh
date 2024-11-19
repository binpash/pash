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

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/max-temp"
INPUTS_DIR="$BENCHMARK_DIR/inputs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="max-temp"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"

INPUTS=(
  # temperatures_200M.txt
  # temperatures_500M.txt
  # temperatures_1G.txt
  temperatures.2015.txt
)

for INPUT in "${INPUTS[@]}"; do
  INPUT_PATH="$INPUTS_DIR/$INPUT"
  S3_URI=$S3_BUCKET_PREFIX/$S3_INPUTS_DIR/$INPUT
  echo "Uploading $INPUT_PATH to $S3_URI"
  aws s3 cp "$INPUT_PATH" "$S3_URI"
done
