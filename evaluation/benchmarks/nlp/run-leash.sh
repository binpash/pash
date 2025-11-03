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
    6_4.sh
    # 6_5.sh
    # 6_7.sh
    # 7_1.sh
    # 7_2.sh
    # 8.2_1.sh
    # 6_1.sh
    # 8_1.sh
    # 8.3_2.sh
    # 8.3_3.sh
  )

for script in "${SCRIPTS[@]}"
do
  echo "ENTRIES=1000 IN=$INFILE OUT=nlp/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}- $PASH_TOP/pa.sh -w${WIDTH} --serverless_exec --parallel_pipelines scripts/$script" 
  { time ENTRIES=1000 IN=$INFILE OUT=nlp/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}- $PASH_TOP/pa.sh --serverless_exec -w${WIDTH} --parallel_pipelines scripts/$script; } 2>outputs/leash-${script}-pash-w${WIDTH}-time.log

  sleep 20
  logs_dir="logs/$script"
  if [ -d "$logs_dir" ]; then
      echo "Removing existing logs directory: $logs_dir"
      rm -rf "$logs_dir"
  fi
  python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir" >outputs/leash-${script}-pash-w${WIDTH}-log-analysis.log
done
