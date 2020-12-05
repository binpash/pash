#!/bin/bash

input_vars_file=$1
output_vars_file=$2

# >&2 echo $input_vars_file

source pash_source_declare_vars.sh $input_vars_file

>&2 echo "Executing script:"
>&2 cat "${@:3}"   
source "${@:3}"
pash_exec_status=$?
declare -p > $output_vars_file
(exit "$pash_exec_status")
