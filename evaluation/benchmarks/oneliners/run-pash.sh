# time IN=oneliners/inputs/3G.txt OUT=oneliners/outputs/nfa-regex-3G-pash-w4-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w4 scripts/nfa-regex.sh

# set -x
# time IN=oneliners/inputs/3G.txt OUT=oneliners/outputs/nfa-regex-3G-pash-w4-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w4 scripts/nfa-regex.sh

# for script in nfa-regex-1.sh nfa-regex-2.sh nfa-regex-3.sh nfa-regex-1-sort.sh nfa-regex-2-sort.sh nfa-regex-3-sort.sh nfa-regex-sort.sh

#!/bin/bash
cd $(dirname $0)

INPUT_SIZE=3G
# for script in sort-sort.sh sort.sh top-n.sh wf.sh set-diff.sh
for script in spell.sh
do
  echo "IN=oneliners/inputs/${INPUT_SIZE}.txt OUT=oneliners/outputs/${script}-${INPUT_SIZE}-pash-w4-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w4 scripts/$script" 
  { time IN=oneliners/inputs/${INPUT_SIZE}.txt OUT=oneliners/outputs/${script}-${INPUT_SIZE}-pash-w4-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w4 scripts/$script; } 2>outputs/${script}-${INPUT_SIZE}-pash-w4-4cpu-time.log
done

# INPUT_SIZE=1G
# for script in bi-grams.sh
# do
#   echo "IN=oneliners/inputs/${INPUT_SIZE}.txt OUT=oneliners/outputs/${script}-${INPUT_SIZE}-pash-w4-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w4 scripts/$script" 
#   { time IN=oneliners/inputs/${INPUT_SIZE}.txt OUT=oneliners/outputs/${script}-${INPUT_SIZE}-pash-w4-4cpu- $PASH_TOP/pa.sh -w4 -d1 scripts/$script; } 2>outputs/${script}-${INPUT_SIZE}-pash-w4-4cpu-time.log
# done