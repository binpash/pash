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

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/nlp"
INPUTS_DIR="$BENCHMARK_DIR/inputs"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"
TIMES_DIR="$BENCHMARK_DIR/times"
SCRIPTS_DIR="$BENCHMARK_DIR/scripts"

S3_BENCHMARK_DIR="nlp"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"
S3_OUTPUTS_DIR="$S3_BENCHMARK_DIR/outputs"

# ENVIRONMENT, MEMORY, SYSTEM, WIDTH
CONFIGS=(
  # Local:16G:Bash:1
  # AWS:2048M:Bash:1
  AWS:2048M:Splash:1
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
    # 1_1.sh:pg-small
    # 2_1.sh:pg-small
    # 2_2.sh:pg-small
    # 3_1.sh:pg-small
    # 3_2.sh:pg-small
    # 3_3.sh:pg-small
    # 4_3b.sh:pg-small
    # 4_3.sh:pg-small
    # 6_1_1.sh:pg-small
    # 6_1_2.sh:pg-small
    # 6_1.sh:pg-small
    # 6_2.sh:pg-small
    # 6_3.sh:pg-small
    # 6_4.sh:pg-small
    # 6_5.sh:pg-small
    # 6_7.sh:pg-small
    # 7_1.sh:pg-small
    # 7_2.sh:pg-small
    # 8_1.sh:pg-small
    # 8.2_1.sh:pg-small
    # 8.2_2.sh:pg-small
    # 8.3_2.sh:pg-small
    # 8.3_3.sh:pg-small
  )
  INPUT_TYPE=".small"
else
  SCRIPTS_INPUTS=(
    # # 1_1_batch.sh:pg
    # 1_1.sh:pg
    # 2_1.sh:pg
    # 2_2.sh:pg
    # 3_1.sh:pg
    # 3_2.sh:pg
    # 3_3.sh:pg
    # 4_3.sh:pg
    # 6_1_1.sh:pg
    # 6_1_2.sh:pg
    # 6_2.sh:pg
    # 6_3.sh:pg
    # 6_4.sh:pg
    # 6_5.sh:pg
    # 6_7.sh:pg
    # 7_1.sh:pg
    # 7_2.sh:pg
    # 8.2_1.sh:pg
    # 6_1.sh:pg
    # 8_1.sh:pg
    8.3_2.sh:pg
  )
  INPUT_TYPE=""
fi

for SCRIPT_INPUT in "${SCRIPTS_INPUTS[@]}"
do
  for CONFIG in "${CONFIGS[@]}"
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
    INPUT_PATH="${INPUTS_DIR}/${INPUT}/"
    OUTPUT_PATH="${OUTPUTS_DIR}/${OUTPUT}"
    TIME_PATH="${TIMES_DIR}/${TIME}"

    S3_INPUT_PATH="${S3_INPUTS_DIR}/${INPUT}/"
    S3_OUTPUT_PATH="${S3_OUTPUTS_DIR}/${OUTPUT}/"

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
      { time IN=$S3_INPUT_PATH OUT=$S3_OUTPUT_PATH ENTRIES=1000 "$PASH_TOP"/pa.sh -w "$WIDTH" -d 1 "$SCRIPT_PATH" --parallel_pipelines --serverless_exec; } 2>"$TIME_PATH"
    fi

    grep real "$TIME_PATH"
  done
done
