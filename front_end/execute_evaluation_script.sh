#!/bin/bash

## We assume that each evaluation script has a sequenital, a
## distributed, and an environment
experiment=$1
directory="../evaluation/"
prefix="${directory}${experiment}"
env_file="${prefix}_env.sh"
seq_script="${prefix}_seq.sh"
distr_script="${prefix}_distr.sh"

echo "Environment:"
cat $env_file
. $env_file
export $(cut -d= -f1 $env_file)

echo "Sequential:"
cat $seq_script
time /bin/bash $seq_script > /tmp/seq_output

echo "Distributed:"
cat $distr_script
time /bin/bash $distr_script
