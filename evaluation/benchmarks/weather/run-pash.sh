#!/bin/bash
cd $(dirname $0)

INFILE=weather/inputs/temperatures.2015.txt
mkdir -p outputs
WIDTH=4

for script in weather-process.sh 
do
  echo "IN=$INFILE OUT=weather/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} scripts/$script" 
  { time IN=$INFILE OUT=weather/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} --parallel_pipelines --parallel_pipelines_limit 4 scripts/$script; } 2>outputs/localecorrect-${script}-pash-w${WIDTH}-4cpu-time.log
done
