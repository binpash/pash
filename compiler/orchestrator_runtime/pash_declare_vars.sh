#!/bin/bash

vars_file="${1?File not given}"

pash_redir_output echo "Writing vars to: $vars_file"

declare -p > "$vars_file"
## KK  2021-11-23 We don't actually need to export functions in the vars file. 
##                We never expand them in the compiler
## declare -f >> "$vars_file"
