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

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/unix50"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="unix50"
S3_OUTPUTS_DIR="$S3_BENCHMARK_DIR/outputs"

CONFIGS=(
  # AWS:2048M:Bash:1
  # AWS:2048M:Splash:1
  # AWS:2048M:Splash:2
  # AWS:2048M:Splash:4
  # AWS:2048M:Splash:8
  # AWS:2048M:Splash:16
  # AWS:2048M:Splash:32
  # AWS:2048M:Splash:64
)


SCRIPTS=(
  # 1.sh to 36.sh autocomplete
  1.sh
  2.sh
  3.sh
  4.sh
  5.sh
  6.sh
  7.sh
  8.sh
  9.sh
  10.sh
  11.sh
  12.sh
  13.sh
  14.sh
  15.sh
  16.sh
  17.sh
  18.sh
  19.sh
  20.sh
  21.sh
  # 22.sh
  23.sh
  24.sh
  25.sh
  26.sh
  # 27.sh
  28.sh
  29.sh
  30.sh
  31.sh
  32.sh
  33.sh
  34.sh
  35.sh
  36.sh
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
        S3_URI=$S3_BUCKET_PREFIX/$S3_OUTPUTS_DIR/$OUTPUT/stdout.txt

        echo "Downloading $S3_URI to $OUTPUT_PATH"

        aws s3 cp  "$S3_URI" "$OUTPUTS_DIR/$OUTPUT"
    done
done
