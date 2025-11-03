#!/bin/bash

cd "$(dirname "$0")" || exit 1

SCRIPT_INPUT_WIDTH=(
    "1.sh:in.csv:16"
    "2.sh:in.csv:16"
    "3.sh:in.csv:16"
    "4.sh:in.csv:16"
)

for SCRIPT_INPUT in "${SCRIPT_INPUT_WIDTH[@]}"; do
    echo "Running benchmark for $SCRIPT_INPUT"
    SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
    INPUT=$(echo "$SCRIPT_INPUT" | cut -d: -f2)
    WIDTH=$(echo "$SCRIPT_INPUT" | cut -d: -f3)
    echo "Running $SCRIPT with input $INPUT and width $WIDTH"
    time IN="covid-mts/inputs/$INPUT" OUT="covid-mts/outputs/$SCRIPT:$INPUT:$WIDTH:hybrid" $PASH_TOP/pa.sh --serverless_exec -w"$WIDTH" scripts/"$SCRIPT"

    sleep 20
    logs_dir="logs/$SCRIPT:$INPUT:$WIDTH"
    if [ -d "$logs_dir" ]; then
        echo "Removing existing logs directory: $logs_dir"
        rm -rf "$logs_dir"
    fi
    python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir"
done
