#!/bin/bash

cd "$(dirname "$0")" || exit 1

SCRIPT=max-temp-process.sh
INPUT=temperatures.2015.txt
# INPUT=temperatures_1G.txt
WIDTH=16

echo  "Running script: $SCRIPT with input: $INPUT and width: $WIDTH"
time IN="max-temp/inputs/$INPUT" OUT="max-temp/outputs/$SCRIPT:$INPUT:$WIDTH:hybrid:" $PASH_TOP/pa.sh --serverless_exec -w"$WIDTH" --parallel_pipelines --parallel_pipelines_limit 16 scripts/"$SCRIPT" 

sleep 10
logs_dir="logs/$SCRIPT:$INPUT:$WIDTH"
if [ -d "$logs_dir" ]; then
    echo "Removing existing logs directory: $logs_dir"
    rm -rf "$logs_dir"
fi
python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir"
