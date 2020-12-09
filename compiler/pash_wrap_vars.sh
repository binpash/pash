#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

input_vars_file=$1
output_vars_file=$2

# pash_redir_output echo $input_vars_file

source "$RUNTIME_DIR/pash_source_declare_vars.sh" $input_vars_file

pash_redir_output echo "Executing script in ${@:3}:"
pash_redir_output cat "${@:3}"
(exit "$pash_previous_exit_status")
source "${@:3}"
pash_exec_status=$?
source "$RUNTIME_DIR/pash_declare_vars.sh" $output_vars_file
(exit "$pash_exec_status")
