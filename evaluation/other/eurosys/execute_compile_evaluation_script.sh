#!/bin/bash

## Necessary to set PASH_TOP
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

execute_seq_flag=0
eager_flag=0
no_task_par_eager_flag=0
auto_split_flag=0
assert_compiler_success=""

while getopts 'senpac' opt; do
    case $opt in
        s) execute_seq_flag=1 ;;
        e) eager_flag=1 ;;
        n) no_task_par_eager_flag=1 ;;
        a) auto_split_flag=1 ;;
        c) assert_compiler_success="--assert_compiler_success" ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

## We assume that each evaluation script has a sequential, a
## distributed, and an environment
microbenchmark=$1
n_in=$2
results_subdir=$3
intermediary_prefix=$4

experiment="${microbenchmark}_${n_in}"

eval_directory="$PASH_TOP/evaluation/"
directory="${eval_directory}/${4}intermediary/"
results="${eval_directory}results/${results_subdir}/"
prefix="${directory}${experiment}"
env_file="${prefix}_env.sh"
funs_file="${prefix}_funs.sh"
seq_script="${prefix}_seq.sh"
input_file="${prefix}.in"

seq_output="${directory}/${microbenchmark}_seq_output"
pash_output="${directory}/${microbenchmark}_pash_output"

echo "Environment:"
# cat $env_file
. $env_file
vars_to_export=$(cut -d= -f1 $env_file)
if [ ! -z "$vars_to_export" ]; then
    export $vars_to_export
fi

## Export necessary functions
if [ -f "$funs_file" ]; then
    source $funs_file
fi

## Redirect the input if there is an input file
stdin_redir="/dev/null"
if [ -f "$input_file" ]; then
    stdin_redir="$(cat "$input_file")"
    echo "Has input file: $stdin_redir"
fi

## TODO: Extend this script to give input to some arguments from stdin.

if [ "$execute_seq_flag" -eq 1 ]; then
    echo "Sequential:"
    cat $seq_script
    cat $stdin_redir | { time /bin/bash $seq_script > $seq_output ; } 2> >(tee "${results}${experiment}_seq.time" >&2)
else
    echo "Not executing sequential..."
fi

## Save the configuration to restore it afterwards
auto_split_opt="--width 1"
config_path_opt=""

if [ "$auto_split_flag" -eq 1 ]; then
    echo "Distributed with auto-split:"
    eager_opt=""
    auto_split_opt="--width ${n_in}"
    distr_result_filename="${results}${experiment}_distr_auto_split.time"
elif [ "$eager_flag" -eq 1 ]; then
    echo "Distributed:"
    eager_opt=""
    distr_result_filename="${results}${experiment}_distr.time"
elif [ "$no_task_par_eager_flag" -eq 1 ]; then
    echo "Distributed with naive (no-task-par) eager:"
    eager_opt=""
    distr_result_filename="${results}${experiment}_distr_no_task_par_eager.time"

    ## Change the configuration
    config_path="/tmp/new-config.yaml"
    config_path_opt="--config_path ${config_path}"
    cat "$PASH_TOP/compiler/config.yaml" > ${config_path}
    sed -i 's/runtime\/eager.sh/runtime\/eager-no-task-par.sh/g' "${config_path}"
else
    echo "Distributed without eager:"
    eager_opt="--no_eager"
    distr_result_filename="${results}${experiment}_distr_no_eager.time"
fi

cat $stdin_redir | { time python3 $PASH_TOP/compiler/pash.py -d 1 --speculation no_spec $assert_compiler_success $eager_opt $auto_split_opt $config_path_opt --output_time $seq_script ; } 1> $pash_output 2> >(tee "${distr_result_filename}" >&2) &&
echo "Checking for equivalence..." &&
diff -s $seq_output $pash_output | head | tee -a "${distr_result_filename}" >&2
