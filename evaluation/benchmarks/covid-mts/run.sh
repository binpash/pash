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

[ -z "$AWS_QUEUE" ] && {
  echo "AWS_QUEUE not set"
  exit
}

mkdir -p {outputs,times}/

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/covid-mts"
INPUTS_DIR="$BENCHMARK_DIR/inputs"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"
TIMES_DIR="$BENCHMARK_DIR/times"
SCRIPTS_DIR="$BENCHMARK_DIR/scripts"

S3_BENCHMARK_DIR="covid-mts"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"
S3_OUTPUTS_DIR="$S3_BENCHMARK_DIR/outputs"

# ENVIRONMENT, MEMORY, SYSTEM, WIDTH
CONFIGS=(
  Local:16G:Bash:1
  # AWS:2048M:Bash:1
  # AWS:2048M:Splash:1
  # AWS:2048M:Splash:2
  # AWS:2048M:Splash:4
  # AWS:2048M:Splash:8
  # AWS:2048M:Splash:16
  # AWS:2048M:Splash:32
  # AWS:2048M:Splash:64
)

if [[ "$*" == *"--small"* ]]
then
  SCRIPTS_INPUTS=(
    1.sh:in_small.csv
    2.sh:in_small.csv
    3.sh:in_small.csv
    4.sh:in_small.csv
    # 5.sh:in_small.csv
  )
  INPUT_TYPE=".small"
else
  SCRIPTS_INPUTS=(
    1.sh:in.csv
    2.sh:in.csv
    3.sh:in.csv
    4.sh:in.csv
    # 5.sh:in.csv
  )
  INPUT_TYPE=""
fi

for CONFIG in "${CONFIGS[@]}"
do
  for SCRIPT_INPUT in "${SCRIPTS_INPUTS[@]}"
  do
    ENVIRONMENT=$(echo "$CONFIG" | cut -d: -f1)
    MEMORY=$(echo "$CONFIG" | cut -d: -f2)
    SYSTEM=$(echo "$CONFIG" | cut -d: -f3)
    WIDTH=$(echo "$CONFIG" | cut -d: -f4)

    SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
    INPUT=$(echo "$SCRIPT_INPUT" | cut -d: -f2)

    OUTPUT="${SCRIPT}__env${ENVIRONMENT}__mem${MEMORY}__sys${SYSTEM}__w${WIDTH}${INPUT_TYPE}.out"
    TIME="${SCRIPT}__env${ENVIRONMENT}__mem${MEMORY}__sys${SYSTEM}__w${WIDTH}${INPUT_TYPE}.time"

    SCRIPT_PATH="${SCRIPTS_DIR}/${SCRIPT}"
    INPUT_PATH="${INPUTS_DIR}/${INPUT}"
    OUTPUT_PATH="${OUTPUTS_DIR}/${OUTPUT}"
    TIME_PATH="${TIMES_DIR}/${TIME}"

    S3_INPUT_PATH="${S3_INPUTS_DIR}/${INPUT}"
    S3_OUTPUT_PATH="${S3_OUTPUTS_DIR}/${OUTPUT}"

    echo "CONFIG: $CONFIG, SCRIPT_INPUT: $SCRIPT_INPUT"

    if [[ $ENVIRONMENT == "Local" ]]
    then
      { time IN=$INPUT_PATH bash "$SCRIPT_PATH" >"$OUTPUT_PATH" ; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Bash" ]]
    then
      echo "Not implemented yet"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Splash" ]]
    then
      echo "Running $SCRIPT with $ENVIRONMENT:$SYSTEM on $S3_INPUT_PATH"
      python3 "$PASH_TOP"/scripts/serverless/delete-log-streams.py
      time IN=$S3_INPUT_PATH "$PASH_TOP"/pa.sh -w "$WIDTH" --serverless_exec "$SCRIPT_PATH" --sls_output "$S3_OUTPUT_PATH"
      { time IN=$S3_INPUT_PATH "$PASH_TOP"/pa.sh -w "$WIDTH" --serverless_exec "$SCRIPT_PATH" --sls_output "$S3_OUTPUT_PATH" ; } 2>"$TIME_PATH"
    fi

    grep real "$TIME_PATH"
  done
done
