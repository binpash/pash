#!/bin/bash

cd $(dirname "$0")

JOB_ID="test-ft-nfa-regex-faulty"
bash test_stateless_faulty.sh 2>&1 >ec2.log

sleep 12

mkdir -p debug
python3 $PASH_TOP/scripts/serverless/utils.py debug