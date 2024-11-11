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

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/unix50"
INPUTS_DIR="$BENCHMARK_DIR/inputs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="unix50"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"

INPUTS=(
  # 1_1M.txt
  # 10_1M.txt
  # 11_1M.txt
  # 12_1M.txt
  # 2_1M.txt
  # 3_1M.txt
  # 4_1M.txt
  # 5_1M.txt
  # 6_1M.txt
  # 7_1M.txt
  # 8_1M.txt
  # 9.1_1M.txt
  # 9.2_1M.txt
  # 9.3_1M.txt
  # 9.4_1M.txt
  # 9.5_1M.txt
  # 9.6_1M.txt
  # 9.7_1M.txt
  # 9.8_1M.txt
  # 9.9_1M.txt

  # 1_1G.txt
  # 10_1G.txt
  # 11_1G.txt
  # 12_1G.txt
  # 2_1G.txt
  # 3_1G.txt
  # 4_1G.txt
  # 5_1G.txt
  # 6_1G.txt
  # 7_1G.txt
  # 8_1G.txt
  # 9.1_1G.txt
  # 9.2_1G.txt
  # 9.3_1G.txt
  # 9.4_1G.txt
  # 9.5_1G.txt
  # 9.6_1G.txt
  # 9.7_1G.txt
  # 9.8_1G.txt
  # 9.9_1G.txt
)

for INPUT in "${INPUTS[@]}"; do
  INPUT_PATH="$INPUTS_DIR/$INPUT"
  S3_URI=$S3_BUCKET_PREFIX/$S3_INPUTS_DIR/$INPUT
  echo "Uploading $INPUT_PATH to $S3_URI"
  aws s3 cp "$INPUT_PATH" "$S3_URI"
done
