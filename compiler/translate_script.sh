#!/bin/bash

seq_script=$1
distr_script=$2

## The second argument is the directory that the parallel script is
## going to output
distr_output_dir=$3

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

## TODO: This currently assumes that there is only one IR in the script
python3.8 $DISH_TOP/compiler/dish.py --output_optimized $seq_script $distr_script &&
cat $distr_script
