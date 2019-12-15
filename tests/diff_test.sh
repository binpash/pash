#!/bin/bash

## We assume that each test script has a sequential version and an
## environment
experiment=$1
distr_output_dir=$2
fan_out=$3
batch_size=$4

directory="${DISH_TOP}/tests/test_scripts/"
prefix="${directory}${experiment}"
env_file="${prefix}_env.sh"
seq_script="${prefix}_seq.sh"
distr_script="/tmp/${experiment}_distr.sh"

seq_output="/tmp/seq_output"
distr_output="/tmp/distr_output_combined"

echo "Environment:"
# cat $env_file
. $env_file
export $(cut -d= -f1 $env_file)

echo "Sequential:"
cat $seq_script | grep -v "rm -f" | grep -v "mkfifo"
/bin/bash $seq_script > $seq_output

echo "Distributed:"
$DISH_TOP/compiler/translate_script.sh $seq_script $distr_script $distr_output_dir $fan_out $batch_size
cat $distr_output_dir/* > $distr_output

echo "Computing output difference... This will take some time..."
## TODO: Solve the discrepancy in the sequential and distributed output
diff -q $seq_output $distr_output
if [[ $? == "0" ]]
then
    echo "Output files are the same :)"
else
    echo "ERROR: Output files are not the same"
    exit 1
fi
