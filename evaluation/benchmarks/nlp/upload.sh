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

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/nlp"
# INPUTS_DIR="$BENCHMARK_DIR/inputs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="nlp"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"

INPUTS_DIRS=(
  pg
  # pg-small
)

# upload small inputs
for INPUT_DIR in "${INPUTS_DIRS[@]}"; do
  for INPUT in $(cat $BENCHMARK_DIR/1000-books.txt); do
    INPUT_PATH="$BENCHMARK_DIR/inputs/$INPUT_DIR/$INPUT"
    S3_URI=$S3_BUCKET_PREFIX/$S3_INPUTS_DIR/$INPUT_DIR/$INPUT
    echo "Uploading $INPUT_PATH to $S3_URI"
    aws s3 cp "$INPUT_PATH" "$S3_URI"
  done
done
