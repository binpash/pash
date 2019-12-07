#!/bin/bash

## The first argument is the prefix of the script to translate
script_prefix=$1

## The second argument is the directory that the parallel script is
## going to output
distr_output_dir=$2

if [ -z "$3" ]
  then
    fan_out=5
else
    fan_out=$3
fi

if [ -z "$4" ]
  then
    batch_size=10000
else
    batch_size=$4
fi

json="/tmp/temp_script.json"
seq_script="${script_prefix}_seq.sh"
distr_script="${script_prefix}_distr.sh"

## TODO: This currently assumes that there is only one IR in the script
python3 main.py $seq_script &&
python3 distr_plan.py /tmp/dish_temp_ir_file1 $distr_script $distr_output_dir $fan_out $batch_size
cat $distr_script
