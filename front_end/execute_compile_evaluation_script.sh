#!/bin/bash

## We assume that each evaluation script has a sequential, a
## distributed, and an environment
experiment=$1
distr_output_dir=$2
cat_at_end=$3
fan_out=$4
batch_size=$5


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

# echo "Sequential:"
# cat $seq_script | grep -v "rm -f" | grep -v "mkfifo"
# { time /bin/bash $seq_script > /tmp/seq_output ; } 2> >(tee "${results}${experiment}_seq.time" >&2)

echo "Distributed:"
if [ "$cat_at_end" = true ] ;
then
    echo "Execute with cat"
    { time { ./translate_script.sh $prefix $distr_output_dir $fan_out $batch_size ; /bin/bash $distr_script ; ./cat_output_files.sh $distr_output_dir > /tmp/distr_output_complete ; } ; } 2> >(tee "${results}${experiment}_compile_cat_distr.time" >&2)
else
    { time { ./translate_script.sh $prefix $distr_output_dir $fan_out $batch_size ; /bin/bash $distr_script ; } ; } 2> >(tee "${results}${experiment}_compile_distr.time" >&2)
fi



