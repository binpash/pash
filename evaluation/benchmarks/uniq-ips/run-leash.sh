#!/bin/bash
cd $(dirname $0)

INFILE=uniq-ips/inputs/logs-popcount-org.txt
mkdir -p outputs
WIDTH=32

for script in uniq-ips.sh
do
  echo "IN=$INFILE OUT=uniq-ips/outputs/${script}-${INPUT_SIZE}-pash-w${WIDTH}-4cpu- $PASH_TOP/pa.sh -w${WIDTH} scripts/$script" 
  { time IN=$INFILE OUT=uniq-ips/outputs/${script}-${INPUT_SIZE}-leash-w${WIDTH}-4cpu- $PASH_TOP/pa.sh -w${WIDTH}  --serverless_exec scripts/$script; } 2>outputs/${script}-leash-w${WIDTH}-4cpu-time.log

  sleep 20
  logs_dir="logs/$SCRIPT:$INPUT:$WIDTH"
  if [ -d "$logs_dir" ]; then
      echo "Removing existing logs directory: $logs_dir"
      rm -rf "$logs_dir"
  fi
  python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir"
done
