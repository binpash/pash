#!/bin/bash

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

BENCHMARK_DIR="${PASH_TOP}/evaluation/benchmarks/covid-mts"
INPUTS_DIR="${BENCHMARK_DIR}/inputs"
OUTPUTS_DIR="${BENCHMARK_DIR}/outputs"
SCRIPTS_DIR="${BENCHMARK_DIR}/scripts"

ENVIRONMENT=Local
MEMORY=16G
SYSTEM=Bash
WIDTH=1
INPUT="${INPUTS_DIR}/in.csv"

SCRIPTS=(
    # "1.sh"
    # "2.sh"
    # "3.sh"
    # "4.sh"
    "5.sh"
)

for SCRIPT in "${SCRIPTS[@]}"
do
    SCRIPT_PATH="${SCRIPTS_DIR}/${SCRIPT}"
    OUTPUT="${OUTPUTS_DIR}/${SCRIPT}__env${ENVIRONMENT}__mem${MEMORY}__sys${SYSTEM}__w${WIDTH}.out"

    echo "Executing $SCRIPT..."

    IN=$INPUT "$SCRIPT_PATH" >"$OUTPUT"
done
