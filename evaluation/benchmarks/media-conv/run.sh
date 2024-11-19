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

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/media-conv"
INPUTS_DIR="$BENCHMARK_DIR/inputs"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"
TIMES_DIR="$BENCHMARK_DIR/times"
SCRIPTS_DIR="$BENCHMARK_DIR/scripts"

S3_BENCHMARK_DIR="media-conv"
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
  )
  INPUT_TYPE=".small"
else
  SCRIPTS_INPUTS=(
    # img_convert.sh:jpg_full/jpg
    to_mp3.sh:wav_full/wav
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
