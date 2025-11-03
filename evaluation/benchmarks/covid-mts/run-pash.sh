#!/bin/bash
cd $(dirname $0)

mkdir -p outputs
WIDTH=2

for script in 1.sh 2.sh 3.sh 4.sh
do
  echo "IN=covid-mts/inputs/in.csv OUT=covid-mts/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} scripts/$script" 
  { time IN=covid-mts/inputs/in.csv OUT=covid-mts/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} scripts/$script -d1 --graphviz pdf; } 2>outputs/localecorrect-${script}-pash-w${WIDTH}-4cpu-time.log
done
