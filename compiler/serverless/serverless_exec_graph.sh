#!/bin/bash
ir_file=$1

response=($(echo "Exec-Graph: $ir_file $declared_functions" | nc -U "$SERVERLESS_SOCKET"))

status=${response[0]} #do something if false
folder_id=${response[1]}
script_id=${response[2]}

if [ "$status" == "ERR" ]; then
    echo "Error in launching serverless execution"
    exit 1
fi