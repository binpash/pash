#!/bin/bash

## We assume that each evaluation script has a sequential, a
## distributed, and an environment
experiment=$1

directory="../evaluation/"
results="${directory}results/"
prefix="${directory}${experiment}"
env_file="${prefix}_env.sh"
seq_script="${prefix}_seq.sh"
distr_script="${prefix}_distr.sh"

echo "Environment:"
# cat $env_file
. $env_file
export $(cut -d= -f1 $env_file)

echo "Sequential:"
cat $seq_script | grep -v "rm -f" | grep -v "mkfifo"
{ time /bin/bash $seq_script > /tmp/seq_output ; } 2> >(tee "${results}${experiment}_seq.time" >&2)

echo "Distributed:"
{ time ./translate_script.sh $seq_script $distr_script ; } 2> >(tee "${results}${experiment}_compile_distr.time" >&2)

## Obsolete as dish.py executes output script
# { time /bin/bash $distr_script ; } 2> >(tee "${results}${experiment}_distr.time" >&2)
## TODO: Move in impl.py and add a flag in main that enables/disables
## concatenating the output files
# { time ./cat_output_files.sh $distr_output_dir > /tmp/distr_output_complete ; } 2> >(tee "${results}${experiment}_cat_distr.time" >&2)




