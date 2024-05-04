#!/bin/bash

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

BENCHMARK_DIR=$PASH_TOP/evaluation/benchmarks/oneliners
INPUTS_DIR=$BENCHMARK_DIR/inputs
OUTPUTS_DIR=$BENCHMARK_DIR/outputs
SCRIPTS_DIR=$BENCHMARK_DIR/scripts

SCRIPTS_INPUTS=(
    # nfa-regex.sh:1M.txt
    # shortest-scripts.sh:all_cmds.txt
    # sort-sort.sh:1M.txt
    # sort.sh:1M.txt
    # spell.sh:1M.txt
    # top-n.sh:1M.txt
    wf.sh:1M.txt
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
