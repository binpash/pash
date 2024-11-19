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
  # AWS:2048M:Lambda_BASH:1
  # AWS:2048M:Splash:1
  # AWS:2048M:Splash:2
  AWS:2048M:Splash:4
  # AWS:2048M:Splash:8
  # AWS:2048M:Splash:16
  # AWS:2048M:Splash:64
  # AWS:2048M:Splash:128
  # AWS:2048M:EC2_BASH:1
  # AWS:2048M:Hybrid:1
  # AWS:2048M:Hybrid:2
  # AWS:2048M:Hybrid:4
  # AWS:2048M:Hybrid:8
  # AWS:2048M:Hybrid:16
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
    # "1.sh:1_20G.txt"
    # "3.sh:1_20G.txt" # Error
    # "2.sh:1_20G.txt" # 4096 Mem usage
    # "4.sh:1_20G.txt" # 4096 Mem usage
    # "5.sh:2_20G.txt"
    # "6.sh:3_20G.txt"
    # "7.sh:4_20G.txt"
    # "8.sh:4_20G.txt"
    # "9.sh:4_20G.txt"
    # "10.sh:4_20G.txt" # sort

    # "11.sh:4_20G.txt" # sort
    12.sh:4_20G.txt
    # "13.sh:5_20G.txt"
    # "14.sh:6_20G.txt" # sort
    # "15.sh:7_20G.txt" 
    # "17.sh:7_20G.txt" # sort
    # "18.sh:8_20G.txt"
    # "19.sh:8_20G.txt"
    # # "20.sh:8_20G.txt"
    # "21.sh:8_20G.txt"
    # "23.sh:9.1_20G.txt"
    # "24.sh:9.2_20G.txt"
    # "25.sh:9.3_20G.txt"
    # "26.sh:9.4_20G.txt"
    # "28.sh:9.6_20G.txt"
    # "29.sh:9.7_20G.txt"
    # "30.sh:9.8_20G.txt"
    # "31.sh:9.9_20G.txt"
    # "32.sh:10_20G.txt"
    # "33.sh:10_20G.txt"
    # "35.sh:11_20G.txt"

    # "34.sh:10_20G.txt"
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
    LOG="${SCRIPT}__env${ENVIRONMENT}__mem${MEMORY}__sys${SYSTEM}__w${WIDTH}${INPUT_TYPE}.log"

    SCRIPT_PATH="${SCRIPTS_DIR}/${SCRIPT}"
    INPUT_PATH="${INPUTS_DIR}/${INPUT}"
    OUTPUT_PATH="${OUTPUTS_DIR}/${OUTPUT}"
    TIME_PATH="${TIMES_DIR}/${TIME}"
    LOG_PATH="${TIMES_DIR}/${LOG}/"

    S3_INPUT_PATH="${S3_INPUTS_DIR}/${INPUT}"
    S3_OUTPUT_PATH="${S3_OUTPUTS_DIR}/${OUTPUT}/" # make sure this is a directory

    echo "CONFIG: $CONFIG, SCRIPT_INPUT: $SCRIPT_INPUT"

    if [[ $ENVIRONMENT == "Local" && $SYSTEM == "Bash" ]]
    then
      { time IN=$INPUT_PATH bash "$SCRIPT_PATH" >"$OUTPUT_PATH" ; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "Local" && $SYSTEM == "PaSh" ]]
    then
      { time IN=$INPUT_PATH "$PASH_TOP"/pa.sh -w "$WIDTH" "$SCRIPT_PATH" >"$OUTPUT_PATH" ; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Bash" ]]
    then
      echo "Not implemented yet"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Lambda_BASH" ]]
    then
      python3 "$PASH_TOP"/scripts/serverless/delete-log-streams.py
      { time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" python3 "$PASH_TOP"/scripts/serverless/run-aws-bash.py "$SCRIPT_PATH-bash" "lambda"; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "EC2_BASH" ]]
    then
      python3 "$PASH_TOP"/scripts/serverless/delete-log-streams.py
      { time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" python3 "$PASH_TOP"/scripts/serverless/run-aws-bash.py "$SCRIPT_PATH-bash" "hybrid"; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Hybrid" ]]
    then
      python3 "$PASH_TOP"/scripts/serverless/delete-log-streams.py
      { time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" "$PASH_TOP"/pa.sh  -d 1 -w "$WIDTH" "$SCRIPT_PATH" --serverless_exec --sls_instance hybrid ; } 2>"$TIME_PATH"
    elif [[ $ENVIRONMENT == "AWS" && $SYSTEM == "Splash" ]]
    then
      python3 "$PASH_TOP"/scripts/serverless/delete-log-streams.py
      { timeout 920 bash -c "time IN=$S3_INPUT_PATH OUT=$"$S3_OUTPUT_PATH" DICT="$S3_INPUTS_DIR/dict.txt" "$PASH_TOP"/pa.sh --graphviz pdf -d 1 -w "$WIDTH" "$SCRIPT_PATH" --serverless_exec" ; } 2>"$TIME_PATH"
      sleep 60
      # sleep 480 # wait for the ports to be recycled
      rm -rf "$LOG_PATH"
      LOG_FOLDER="$LOG_PATH" python3 "$PASH_TOP"/scripts/serverless/get-log-streams.py
    fi

    grep real "$TIME_PATH"
  done
done
