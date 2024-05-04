#!/bin/bash

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

BENCHMARK_DIR=$PASH_TOP/evaluation/benchmarks/unix50
INPUTS_DIR=$BENCHMARK_DIR/inputs
OUTPUTS_DIR=$BENCHMARK_DIR/outputs
SCRIPTS_DIR=$BENCHMARK_DIR/scripts

# We don't run {22,27}.sh because they are empty
SCRIPTS_INPUTS=(
    # "01.sh:01.txt"
    # "02.sh:01.txt"
    # "03.sh:01.txt"
    # "04.sh:01.txt"
    # "05.sh:02.txt"
    # "06.sh:03.txt"
    # "07.sh:04.txt"
    # "08.sh:04.txt"
    # "09.sh:04.txt"
    # "10.sh:04.txt"
    # "11.sh:04.txt"
    # "12.sh:04.txt"
    # "13.sh:05.txt"
    # "14.sh:06.txt"
    # "15.sh:07.txt"
    # "16.sh:07.txt"
    # "17.sh:07.txt"
    # "18.sh:08.txt"
    # "19.sh:08.txt"
    # "20.sh:08.txt"
    # "21.sh:08.txt"
    # "23.sh:09.1.txt"
    # "24.sh:09.2.txt"
    # "25.sh:09.3.txt"
    # "26.sh:09.4.txt"
    # "28.sh:09.6.txt"
    # "29.sh:09.7.txt"
    # "30.sh:09.8.txt"
    # "31.sh:09.9.txt"
    # "32.sh:10.txt"
    # "33.sh:10.txt"
    # "34.sh:10.txt"
    # "35.sh:11.txt"
    # "36.sh:11.txt"
)

ENVIRONMENT=Local
MEMORY=16G
SYSTEM=Bash
WIDTH=1

for SCRIPT_INPUT in "${SCRIPTS_INPUTS[@]}"
do
  SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
  SCRIPT_PATH="${SCRIPTS_DIR}/${SCRIPT}"
  INPUT=${INPUTS_DIR}/$(echo "$SCRIPT_INPUT" | cut -d: -f2)
  OUTPUT="${OUTPUTS_DIR}/${SCRIPT}__env${ENVIRONMENT}__mem${MEMORY}__sys${SYSTEM}__w${WIDTH}.out"

  IN=$INPUT "$SCRIPT_PATH" >"$OUTPUT"
done
