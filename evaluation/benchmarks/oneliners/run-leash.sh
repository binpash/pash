#!/bin/bash

cd "$(dirname "$0")" || exit 1

SCRIPT_INPUT_WIDTH=(
    "sort.sh:1G.txt:16"
    "sort-sort.sh:3G.txt:16"
    "wf.sh:3G.txt:16"
    "top-n.sh:3G.txt:16"
    "set-diff.sh:3G.txt:16"
    "bi-grams.sh:1G.txt:16"
    "spell.sh:3G.txt:16"
    "nfa-regex.sh:3G.txt:128"
)

for SCRIPT_INPUT in "${SCRIPT_INPUT_WIDTH[@]}"; do
    echo "Running benchmark for $SCRIPT_INPUT"
    SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
    INPUT=$(echo "$SCRIPT_INPUT" | cut -d: -f2)
    WIDTH=$(echo "$SCRIPT_INPUT" | cut -d: -f3)
    echo "Running $SCRIPT with input $INPUT and width $WIDTH"
    time IN="oneliners/inputs/$INPUT" OUT="oneliners/outputs/$SCRIPT:$INPUT:$WIDTH:leash" DICT="oneliners/inputs/dict.txt" $PASH_TOP/pa.sh --serverless_exec -w"$WIDTH" scripts/"$SCRIPT"

    sleep 20
    logs_dir="logs/$SCRIPT:$INPUT:$WIDTH"
    if [ -d "$logs_dir" ]; then
        echo "Removing existing logs directory: $logs_dir"
        rm -rf "$logs_dir"
    fi
    python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir"
done