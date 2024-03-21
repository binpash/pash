#!/bin/bash

vars_file="${1?File not given}"

# pash_redir_output echo "Writing vars to: $vars_file"

echo "cd ${PWD}" > "$vars_file"
declare -p >> "$vars_file"
declare -f >> "$vars_file"
