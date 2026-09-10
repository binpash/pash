#!/bin/bash
cd $(dirname $0)

INFILE=agent/inputs/logs/
mkdir -p outputs
WIDTH=1
N_CPU=$(nproc)

SCRIPTS=(
  terminal-bench-log-summary.sh
)

for script in "${SCRIPTS[@]}"
do
  echo "IN=$INFILE OUT=agent/outputs/${script}-pash-w${WIDTH}-${N_CPU}cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} --parallel_pipelines --parallel_pipelines_limit $N_CPU scripts/$script"
  { time IN=$INFILE OUT=agent/outputs/${script}-pash-w${WIDTH}-${N_CPU}cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} --parallel_pipelines --parallel_pipelines_limit "$N_CPU" scripts/$script; } &>outputs/localecorrect-${script}-pash-w${WIDTH}-${N_CPU}cpu-time.log
done
