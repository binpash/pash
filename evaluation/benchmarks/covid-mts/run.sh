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

mkdir -p {fifos,outputs,times}/

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/covid-mts"
INPUTS_DIR="$BENCHMARK_DIR/inputs"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"
TIMES_DIR="$BENCHMARK_DIR/times"
SCRIPTS_DIR="$BENCHMARK_DIR/scripts"

S3_BENCHMARK_DIR="covid-mts"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"
S3_OUTPUTS_DIR="$S3_BENCHMARK_DIR/outputs"

SERVERLESS_SCRIPTS_DIR="$PASH_TOP/scripts/serverless"
SERVERLESS_RUNTIME_DIR="$PASH_TOP/runtime/serverless/aws"

FIFOS_DIR="$BENCHMARK_DIR/fifos"

# ENVIRONMENT, MEMORY, SYSTEM, WIDTH
CONFIGS=(
  Local:16G:Bash:1
  # AWS:2048M:Bash:1
  # AWS:2048M:Splash:1
  # AWS:2048M:Splash:2
  # AWS:2048M:Splash:4
  # AWS:2048M:Splash:8
  # AWS:2048M:Splash:16
  # AWS:2048M:Hybrid:1
  # AWS:2048M:Hybrid:2
  # AWS:2048M:Hybrid:4
  # AWS:2048M:Hybrid:8
  # AWS:2048M:Hybrid:16
  # AWS:2048M:EC2_BASH:1
  # AWS:2048M:Lambda_BASH:1
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
    1.sh:in_1G.csv
    2.sh:in_1G.csv
    3.sh:in_1G.csv
    4.sh:in_1G.csv
    # 5.sh:in_1G.csv
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
    INPUT_PATH="${INPUTS_DIR}/${INPUT}"
    OUTPUT_PATH="${OUTPUTS_DIR}/${OUTPUT}"
    TIME_PATH="${TIMES_DIR}/${TIME}"

    S3_INPUT_PATH="${S3_INPUTS_DIR}/${INPUT}"
    S3_OUTPUT_PATH="${S3_OUTPUTS_DIR}/${OUTPUT}" # make sure this is a directory

    echo "CONFIG: $CONFIG, SCRIPT_INPUT: $SCRIPT_INPUT"

    if [[ $ENVIRONMENT == "Local" && $SYSTEM == "Bash" ]]
    then
      IN_FIFO_PATH="$FIFOS_DIR"/infifo
      OUT_FIFO_PATH="$FIFOS_DIR"/outfifo
      rm -f "$FIFOS_DIR"/*
      mkfifo "$IN_FIFO_PATH" "${OUT_FIFO_PATH}stdout.txt" # Temp fix, to match with name used in scripts
      python3 "$SERVERLESS_RUNTIME_DIR/s3-get-object.py" "$S3_INPUT_PATH" "$IN_FIFO_PATH" &
      { time IN="$IN_FIFO_PATH" OUT="$OUT_FIFO_PATH" bash "$SCRIPT_PATH" ; } 2>"$TIME_PATH" &
      python3 "$SERVERLESS_RUNTIME_DIR/s3-put-object.py" "$S3_OUTPUT_PATH" "${OUT_FIFO_PATH}stdout.txt" dummy # dummy object_id, not used for bash
    elif [[ $ENVIRONMENT == "Local" && $SYSTEM == "PaSh" ]]
    then
      { time IN=$INPUT_PATH "$PASH_TOP"/pa.sh -w "$WIDTH" "$SCRIPT_PATH" >"$OUTPUT_PATH" ; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Bash" ]]
    then
      echo "Not implemented yet"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Lambda_BASH" ]]
    then
      python3 "$SERVERLESS_SCRIPTS_DIR/delete-log-streams.py"
      { time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" python3 "$SERVERLESS_SCRIPTS_DIR/run-aws-bash.py" "$SCRIPT_PATH-bash" "lambda"; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "EC2_BASH" ]]
    then
      python3 "$SERVERLESS_SCRIPTS_DIR/delete-log-streams.py"
      { time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" python3 "$SERVERLESS_SCRIPTS_DIR/run-aws-bash.py" "$SCRIPT_PATH-bash" "hybrid"; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Hybrid" ]]
    then
      python3 "$SERVERLESS_SCRIPTS_DIR/delete-log-streams.py"
      { time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" "$PASH_TOP"/pa.sh  -d 1 -w "$WIDTH" "$SCRIPT_PATH" --serverless_exec --sls_instance hybrid ; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Splash" ]]
    then
      python3 "$SERVERLESS_SCRIPTS_DIR/delete-log-streams.py"
      { time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" "$PASH_TOP"/pa.sh  -d 1 -w "$WIDTH" "$SCRIPT_PATH" --serverless_exec ; } 2>"$TIME_PATH"
    fi

    grep real "$TIME_PATH"
  done
done
