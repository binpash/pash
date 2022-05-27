#!/bin/bash


microbenchmark=$1
n_in=$2
results_subdir="gnu_parallel"

experiment="${microbenchmark}_${n_in}"

DISH_TOP=${DISH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

eval_directory="../evaluation/"
intermediary_directory="${eval_directory}/intermediary/"
results="${eval_directory}results/${results_subdir}/"
prefix="${intermediary_directory}${experiment}_gnu_parallel"

mkdir -p results

env_file="${prefix}_env.sh"
funs_file="${prefix}_funs.sh"
gnu_parallel_script="${prefix}.sh"

gnu_parallel_scripts_dir="${eval_directory}/gnu_parallel_benchmarks/"
microbenchmarks_dir="${eval_directory}/microbenchmarks/"

## Generate the intermediary gnu parallel scripts
python3 generate_gnu_parallel_intermediary_script.py "${gnu_parallel_scripts_dir}" "${microbenchmarks_dir}" \
        "${microbenchmark}" "${n_in}" "${intermediary_directory}" || 
{ echo 'GNU parallel script generation failed' ; exit 1; }

seq_output="${intermediary_directory}/${microbenchmark}_seq_output"
gnu_parallel_output="${intermediary_directory}/${microbenchmark}_gnu_parallel_output"

echo "Environment:"
cat "$env_file"
. "$env_file"
export "$(cut -d= -f1 "$env_file")"

## Export necessary functions
if [ -f "$funs_file" ]; then
    source "$funs_file"
fi

gnu_parallel_result_filename="${results}${experiment}_gnu_parallel.time"

echo "GNU Parallel:"
cat "$gnu_parallel_script"
{ time /bin/bash "$gnu_parallel_script" > "$gnu_parallel_output" ; } 2> >(tee "$gnu_parallel_result_filename" >&2)

echo "Checking for equivalence..."
diff -s "$seq_output" "$gnu_parallel_output" | tee -a "$gnu_parallel_result_filename"

