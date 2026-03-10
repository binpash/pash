#!/bin/bash
cd $(dirname $0)

INFILE=uniq-ips/inputs/logs-popcount-org-large.txt
mkdir -p outputs
WIDTH=4

for script in uniq-ips.sh-bash
do
  echo "IN=$INFILE OUT=uniq-ips/outputs/${script}-${INPUT_SIZE}-bash- bash scripts/$script" 
  { time IN=$INFILE OUT=uniq-ips/outputs/${script}-${INPUT_SIZE}-bash- bash scripts/$script; } 2>outputs/${script}-bash-time.log
done
