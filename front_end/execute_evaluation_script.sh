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
cat $env_file
. $env_file
export $(cut -d= -f1 $env_file)

echo "Sequential:"
cat $seq_script | grep -v "rm -f" | grep -v "mkfifo"
{ time /bin/bash $seq_script > /tmp/seq_output ; } 2> "${results}${experiment}_seq.time"

echo "Distributed:"
cat $distr_script | grep -v "rm -f" | grep -v "mkfifo"
{ time /bin/bash $distr_script ; } 2> "${results}${experiment}_distr.time"
