#!/bin/bash

input_vars_file=$1
output_vars_file=$2

# >&2 echo $input_vars_file

source pash_source_declare_vars.sh $input_vars_file
source "${@:3}"
declare -p > $output_vars_file