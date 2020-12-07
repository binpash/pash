#!/bin/bash

vars_file="${1?File not given}"

declare -p > "$vars_file"
declare -f >> "$vars_file"