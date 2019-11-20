#!/bin/bash

## The first argument is the prefix of the script to translate
script_prefix=$1

## The second argument is the directory that the parallel script is
## going to output
distr_output_dir=$2

json="/tmp/temp_script.json"
seq_script="${script_prefix}_seq.sh"
distr_script="${script_prefix}_distr.sh"

## First parse the script_to_json
../parser/parse_to_json.native $seq_script > $json

# echo "Translate:"
# cat $seq_script

## TODO: This currently assumes that there is only one IR in the script
python3 main.py $json &&
python3 distr_plan.py /tmp/dish_temp_ir_file1 $distr_script $distr_output_dir
