#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

input_vars_file=$1
output_vars_file=$2

# >&2 echo $input_vars_file

source "$RUNTIME_DIR/pash_source_declare_vars.sh" $input_vars_file

>&2 echo "Executing script in ${@:3}:"
>&2 cat "${@:3}"   
source "${@:3}"
pash_exec_status=$?
source "$RUNTIME_DIR/pash_declare_vars.sh" $output_vars_file
(exit "$pash_exec_status")
