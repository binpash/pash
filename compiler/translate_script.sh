#!/bin/bash

seq_script=$1
distr_script=$2

## TODO: This currently assumes that there is only one IR in the script
python3.8 $DISH_TOP/compiler/dish.py --output_optimized $seq_script $distr_script &&
cat $distr_script
