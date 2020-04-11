#!/bin/bash

execute_seq_flag=0

while getopts 's' opt; do
    case $opt in
        s) execute_seq_flag=1 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

## We assume that each evaluation script has a sequential, a
## distributed, and an environment
experiment=$1
results_subdir=$2

DISH_TOP=${DISH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

eval_directory="../evaluation/"
directory="${eval_directory}/intermediary/"
results="${eval_directory}results/${results_subdir}/"
prefix="${directory}${experiment}"
env_file="${prefix}_env.sh"
funs_file="${prefix}_funs.sh"
seq_script="${prefix}_seq.sh"
distr_script="${prefix}_distr.sh"

echo "Environment:"
# cat $env_file
. $env_file
export $(cut -d= -f1 $env_file)

## Export necessary functions
if [ -f $funs_file ]; then
    source $funs_file
fi

if [ "$execute_seq_flag" -eq 1 ]; then
    echo "Sequential:"
    cat $seq_script
    { time /bin/bash $seq_script > /tmp/seq_output ; } 2> >(tee "${results}${experiment}_seq.time" >&2)
else
    echo "Not executing sequential..."
fi

echo "Distributed:"
{ time python3.8 $DISH_TOP/compiler/dish.py --output_optimized --output_time $seq_script $distr_script ; } 2> >(tee "${results}${experiment}_distr.time" >&2)




