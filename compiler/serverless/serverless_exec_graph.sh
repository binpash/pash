#!/bin/bash
ir_file=$1
# python3 $PASH_TOP/compiler/serverless/serverless_executor.py $1 $declared_functions

response=($(echo "Exec-Graph: $ir_file $declared_functions" | nc -U "$SERVERLESS_SOCKET"))

status=${response[0]} #do something if false
execution_random_id=${response[1]}

if [ "$status" == "ERR" ]; then
    echo "Error in launching serverless execution"
    exit 1
fi

python3 $PASH_TOP/compiler/serverless/serverless_wait_msg.py $execution_random_id