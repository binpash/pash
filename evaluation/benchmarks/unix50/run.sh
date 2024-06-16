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

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/unix50"
INPUTS_DIR="$BENCHMARK_DIR/inputs"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"
TIMES_DIR="$BENCHMARK_DIR/times"
SCRIPTS_DIR="$BENCHMARK_DIR/scripts"

S3_BENCHMARK_DIR="unix50"
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
    # "1.sh:1_1M.txt"
    # # "2.sh:1_1M.txt"
    # "3.sh:1_1M.txt"
    # "4.sh:1_1M.txt"
    # "5.sh:2_1M.txt"
    # "6.sh:3_1M.txt"
    # "7.sh:4_1M.txt"
    # "8.sh:4_1M.txt"
    # "9.sh:4_1M.txt"
    # # "10.sh:4_1M.txt"
    # # "11.sh:4_1M.txt"
    # "13.sh:5_1M.txt"
    # # "14.sh:6_1M.txt"
    # "15.sh:7_1M.txt"
    # "17.sh:7_1M.txt"
    # "18.sh:8_1M.txt"
    # "19.sh:8_1M.txt"
    # "20.sh:8_1M.txt"
    # # "21.sh:8_1M.txt"
    # "23.sh:9.1_1M.txt"
    # "24.sh:9.2_1M.txt"
    # "25.sh:9.3_1M.txt"
    # "26.sh:9.4_1M.txt"
    # "28.sh:9.6_1M.txt"
    # "29.sh:9.7_1M.txt"
    # "30.sh:9.8_1M.txt"
    # "31.sh:9.9_1M.txt"
    # "32.sh:10_1M.txt"
    # "33.sh:10_1M.txt"
    # "35.sh:11_1M.txt"
  )
  INPUT_TYPE=".small"
else
  SCRIPTS_INPUTS=(
    # "1.sh:1_1G.txt"
    # # "2.sh:1_1G.txt"
    # "3.sh:1_1G.txt"
    "4.sh:1_1G.txt"
    # "5.sh:2_1G.txt"
    # "6.sh:3_1G.txt"
    # "7.sh:4_1G.txt"
    # "8.sh:4_1G.txt"
    # "9.sh:4_1G.txt"
    # # "10.sh:4_1G.txt"
    # # "11.sh:4_1G.txt"
    # "13.sh:5_1G.txt"
    # # "14.sh:6_1G.txt"
    # "15.sh:7_1G.txt"
    # "17.sh:7_1G.txt"
    # "18.sh:8_1G.txt"
    # "19.sh:8_1G.txt"
    # "20.sh:8_1G.txt"
    # # "21.sh:8_1G.txt"
    # "23.sh:9.1_1G.txt"
    # "24.sh:9.2_1G.txt"
    # "25.sh:9.3_1G.txt"
    # "26.sh:9.4_1G.txt"
    # "28.sh:9.6_1G.txt"
    # "29.sh:9.7_1G.txt"
    # "30.sh:9.8_1G.txt"
    # "31.sh:9.9_1G.txt"
    # "32.sh:10_1G.txt"
    # "33.sh:10_1G.txt"
    # "35.sh:11_1G.txt"
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
    S3_OUTPUT_PATH="${S3_OUTPUTS_DIR}/${OUTPUT}/"

    echo "CONFIG: $CONFIG, SCRIPT_INPUT: $SCRIPT_INPUT"

    if [[ $ENVIRONMENT == "Local" ]]
    then
      { time IN=$INPUT_PATH bash "$SCRIPT_PATH" $INPUT_PATH >"$OUTPUT_PATH" ; } 2>"$TIME_PATH"
      # IN=$INPUT_PATH bash "$SCRIPT_PATH" $INPUT_PATH ;

    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Bash" ]]
    then
      echo "Not implemented yet"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Splash" ]]
    then
      python3 "$PASH_TOP"/scripts/serverless/delete-log-streams.py
      { time IN=$S3_INPUT_PATH OUT=$S3_OUTPUT_PATH "$PASH_TOP"/pa.sh -w "$WIDTH" -d 1 --serverless_exec "$SCRIPT_PATH" $S3_INPUT_PATH ; } 2>"$TIME_PATH"
    fi

    grep real "$TIME_PATH"
  done
done
