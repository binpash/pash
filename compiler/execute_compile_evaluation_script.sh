#!/bin/bash

execute_seq_flag=0
eager_flag=0
no_task_par_eager_flag=0

while getopts 'sen' opt; do
    case $opt in
        s) execute_seq_flag=1 ;;
        e) eager_flag=1 ;;
        n) no_task_par_eager_flag=1 ;;
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

## Save the configuration to restore it afterwards
cat config.yaml > /tmp/backup-config.yaml

if [ "$eager_flag" -eq 1 ]; then
    echo "Distributed:"
    eager_opt=""
    distr_result_filename="${results}${experiment}_distr.time"
elif [ "$no_task_par_eager_flag" -eq 1 ]; then
    echo "Distributed with naive (no-task-par) eager:"
    eager_opt=""
    distr_result_filename="${results}${experiment}_distr_no_task_par_eager.time"

    ## Change the configuration
    sed -i 's/tools\/eager/tools\/eager-no-task-par.sh/g' config.yaml
else
    echo "Distributed without eager:"
    eager_opt="--no_eager"
    distr_result_filename="${results}${experiment}_distr_no_eager.time"
fi

{ time python3.8 $DISH_TOP/compiler/dish.py --output_optimized $eager_opt --output_time $seq_script $distr_script ; } 2> >(tee "${distr_result_filename}" >&2)

cat /tmp/backup-config.yaml > config.yaml

echo "Checking for equivalence..."
diff -s /tmp/seq_output /tmp/distr_output/0 | tee -a "${distr_result_filename}"

