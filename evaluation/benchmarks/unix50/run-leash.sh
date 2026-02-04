#!/bin/bash

cd "$(dirname "$0")" || exit 1

if [[ "$*" == *"--small"* ]]
then
  SCRIPTS_INPUTS_WIDTH=(
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
  SCRIPTS_INPUTS_WIDTH=(
    "1.sh:1_20G.txt:16"
    "2.sh:1_5G.txt:16" 
    "4.sh:1_5G.txt:16" 
    "5.sh:2_20G.txt:16"
    #"6.sh:3_5G.txt:32" # still fails
    "7.sh:4_20G.txt:16"
    "8.sh:4_20G.txt:16"
    "9.sh:4_20G.txt:16"
    "10.sh:4_20G.txt:16" 
    "11.sh:4_20G.txt:16"
    "13.sh:5_20G.txt:16"
    "14.sh:6_5G.txt:32" 
    "15.sh:7_20G.txt:16" 
    "17.sh:7_20G.txt:16"
    "18.sh:8_20G.txt:16"
    "19.sh:8_20G.txt:16"
    "21.sh:8_5G.txt:32"
    "23.sh:9.1_20G.txt:16"
    "24.sh:9.2_20G.txt:16"
    "25.sh:9.3_20G.txt:16"
    "26.sh:9.4_20G.txt:16"
    "28.sh:9.6_20G.txt:16"
    "29.sh:9.7_20G.txt:16"
    "30.sh:9.8_20G.txt:16"
    "31.sh:9.9_5G.txt:16"
    "32.sh:10_20G.txt:16"
    "33.sh:10_20G.txt:16"
    "35.sh:11_20G.txt:16"
    "36.sh:11_20G.txt:16"
  )
  INPUT_TYPE=""
fi

OUTPUT_DIR="outputs-s3-smart-precheck-complete"
mkdir -p $OUTPUT_DIR

for SCRIPT_INPUT in "${SCRIPTS_INPUTS_WIDTH[@]}"
do
  SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
  INPUT=$(echo "$SCRIPT_INPUT" | cut -d: -f2)
  WIDTH=$(echo "$SCRIPT_INPUT" | cut -d: -f3)

  { time IN=unix50/inputs/$INPUT OUT=unix50/outputs/$SCRIPT:$INPUT:$WIDTH:hybrid "$PASH_TOP"/pa.sh -w "$WIDTH" scripts/"$SCRIPT" --serverless_exec; } 2>$OUTPUT_DIR/leash-$SCRIPT-$INPUT-$WIDTH.time.log
   
  sleep 20
  logs_dir="logs/$SCRIPT:$INPUT:$WIDTH"
  if [ -d "$logs_dir" ]; then
      echo "Removing existing logs directory: $logs_dir"
      rm -rf "$logs_dir"
  fi
  python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir" >$OUTPUT_DIR/leash-$SCRIPT-$INPUT-$WIDTH.analysis.log
done
