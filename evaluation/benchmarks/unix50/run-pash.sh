#!/bin/bash
cd $(dirname $0)

SCRIPTS_INPUTS=(
    # "1.sh:1_20G.txt"
    # "3.sh:1_20G.txt" # Error
    # "2.sh:1_5G.txt" # 4096 Mem usage
    # "4.sh:1_5G.txt" # 4096 Mem usage
    # "5.sh:2_20G.txt"
    # "6.sh:3_20G.txt"
    # "7.sh:4_20G.txt"
    # "8.sh:4_20G.txt"
    # "9.sh:4_20G.txt"
    # "10.sh:4_20G.txt" # sort
    # "11.sh:4_20G.txt" # sort
    # "12.sh:4_20G.txt"
    # "13.sh:5_20G.txt"
    # "14.sh:6_5G.txt" # sort
    # "15.sh:7_20G.txt" 
    # "17.sh:7_20G.txt" # sort
    # "18.sh:8_20G.txt"
    # "19.sh:8_20G.txt"
    # "21.sh:8_5G.txt"
    # "23.sh:9.1_20G.txt"
    # "24.sh:9.2_20G.txt"
    # "25.sh:9.3_20G.txt"
    # "26.sh:9.4_20G.txt"
    # "28.sh:9.6_20G.txt"
    # "29.sh:9.7_20G.txt"
    # "30.sh:9.8_20G.txt"
    # "31.sh:9.9_5G.txt"
    # "32.sh:10_20G.txt"
    # "33.sh:10_20G.txt"
    # "35.sh:11_20G.txt"
    # "34.sh:10_20G.txt"
    "36.sh:11_20G.txt"
  )

mkdir -p outputs
# for script in nfa-regex-1.sh nfa-regex-2.sh nfa-regex-3.sh nfa-regex-1-sort.sh nfa-regex-2-sort.sh nfa-regex-3-sort.sh nfa-regex-sort.sh
for script_input in "${SCRIPTS_INPUTS[@]}"
do
  IFS=':' read -r script INPUT <<< "$script_input"
  echo "IN=unix50/inputs/${INPUT} OUT=unix50/outputs/${script}-${INPUT}-pash-w4-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w4 scripts/$script" 
  { time IN=unix50/inputs/${INPUT} OUT=unix50/outputs/${script}-${INPUT}-pash-w4-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w4 scripts/$script; } 2>outputs/${script}-${INPUT}-pash-w4-4cpu-time.log
done