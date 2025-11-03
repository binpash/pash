#!/bin/bash
cd $(dirname $0)

INFILE=nlp/inputs/pg_heavy/
mkdir -p outputs
WIDTH=1

SCRIPTS=(
    # 1_1.sh
    # 2_1.sh
    # 2_2.sh
    # 3_1.sh
    # 3_2.sh
    # 3_3.sh
    # 4_3.sh
    # 6_1_1.sh
    # 6_1_2.sh
    # 6_2.sh
    # 6_3.sh
    # 6_4.sh
    # 6_5.sh
    # 6_7.sh
    # 7_1.sh
    # 7_2.sh
    # 8.2_1.sh
    # 6_1.sh
    # 8_1.sh
    # 8.3_2.sh
    8.3_3.sh
  )

for script in "${SCRIPTS[@]}"
do
  echo "ENTRIES=1000 IN=$INFILE OUT=nlp/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} --parallel_pipelines --parallel_pipelines_limit 4 scripts/$script" 
  { time ENTRIES=1000 IN=$INFILE OUT=nlp/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} --parallel_pipelines --parallel_pipelines_limit 4 scripts/$script; } 2>outputs/localecorrect-${script}-pash-w${WIDTH}-4cpu-time.log
done
