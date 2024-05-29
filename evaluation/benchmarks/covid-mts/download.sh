#!/bin/bash

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

[ -z "$AWS_ACCOUNT_ID" ] && {
  echo "AWS_ACCOUNT_ID not set"
  exit
}

[ -z "$AWS_BUCKET" ] && {
  echo "AWS_BUCKET not set"
  exit
}

mkdir -p outputs/

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/covid-mts"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="covid-mts"
S3_OUTPUTS_DIR="$S3_BENCHMARK_DIR/outputs"

CONFIGS=(
  # AWS:2048M:Bash:1
  AWS:2048M:Splash:1
  # AWS:2048M:Splash:2
  # AWS:2048M:Splash:4
  # AWS:2048M:Splash:8
  # AWS:2048M:Splash:16
  # AWS:2048M:Splash:32
  # AWS:2048M:Splash:64
)

SCRIPTS=(
  1.sh
  2.sh
  3.sh
  4.sh
  # 5.sh
)

if [[ "$*" == *"--small"* ]]
then
    INPUT_TYPE=".small"
else
    INPUT_TYPE=""
fi

for CONFIG in "${CONFIGS[@]}"
do
    for SCRIPT in "${SCRIPTS[@]}"
    do
        ENVIRONMENT=$(echo "$CONFIG" | cut -d: -f1)
        MEMORY=$(echo "$CONFIG" | cut -d: -f2)
        SYSTEM=$(echo "$CONFIG" | cut -d: -f3)
        WIDTH=$(echo "$CONFIG" | cut -d: -f4)

        OUTPUT="${SCRIPT}__env${ENVIRONMENT}__mem${MEMORY}__sys${SYSTEM}__w${WIDTH}${INPUT_TYPE}.out"

        OUTPUT_PATH="$OUTPUTS_DIR/$OUTPUT"
        S3_URI=$S3_BUCKET_PREFIX/$S3_OUTPUTS_DIR/$OUTPUT

        echo "Downloading $S3_URI to $OUTPUT_PATH"

        aws s3 cp  "$S3_URI" "$OUTPUTS_DIR/$OUTPUT"
    done
done
