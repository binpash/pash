#!/bin/bash

seq_script=$1
distr_script=$2
DISH_TOP=${DISH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
## TODO: This currently assumes that there is only one IR in the script
python3.8 $DISH_TOP/compiler/dish.py --output_optimized $seq_script $distr_script &&
cat $distr_script
