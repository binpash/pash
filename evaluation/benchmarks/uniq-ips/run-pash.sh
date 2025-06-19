#!/bin/bash
cd $(dirname $0)

INFILE=uniq-ips/inputs/logs-popcount-org.txt
mkdir -p outputs
WIDTH=4

for script in uniq-ips.sh
do
  echo "IN=$INFILE OUT=uniq-ips/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} scripts/$script" 
  { time IN=$INFILE OUT=uniq-ips/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w${WIDTH} scripts/$script; } 2>outputs/${script}-pash-w${WIDTH}-4cpu-time.log
done
