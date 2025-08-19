#!/bin/bash

cd "$(dirname "$0")" || exit 1

# Benchmarks
SCRIPT_INPUT_WIDTH=(
    # "sort.sh:3G.txt:16"
    "wf.sh:3G.txt:16"
    "top-n.sh:3G.txt:16"
    # "spell.sh:3G.txt:16"
)

# # Testing
# SCRIPT_INPUT_WIDTH=(
#     "sort.sh:1M.txt:2"
# )

# benchmarks
# time IN=oneliners/inputs/1M.txt OUT=oneliners/outputs/sort.sh.1M.w64 $PASH_TOP/pa.sh --serverless_exec -d1 -w2 scripts/sort.sh
# time OUT=oneliners/outputs/nfa-regex-3G-leash-w128-lambda-only- $PASH_TOP/pa.sh --serverless_exec -w128 scripts/nfa-regex.sh

# W=16
# time IN=oneliners/inputs/3G.txt OUT=oneliners/outputs/sort.sh:$W:hybrid $PASH_TOP/pa.sh --serverless_exec -w$W scripts/sort.sh

for SCRIPT_INPUT in "${SCRIPT_INPUT_WIDTH[@]}"; do
    echo ">>> Running benchmark for $SCRIPT_INPUT"
    SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
    INPUT=$(echo "$SCRIPT_INPUT" | cut -d: -f2)
    WIDTH=$(echo "$SCRIPT_INPUT" | cut -d: -f3)
    echo "Running $SCRIPT with input $INPUT and width $WIDTH"
    time IN="oneliners/inputs/$INPUT" DICT="oneliners/inputs/dict.txt" OUT="oneliners/outputs/$SCRIPT:$INPUT:$WIDTH:hybrid" $PASH_TOP/pa.sh --serverless_exec -w"$WIDTH" scripts/"$SCRIPT"

    sleep 10
    logs_dir="logs/$SCRIPT:$INPUT:$WIDTH"
    if [ -d "$logs_dir" ]; then
        echo "Removing existing logs directory: $logs_dir"
        rm -rf "$logs_dir"
    fi
    python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir"
done
