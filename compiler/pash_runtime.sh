#!/bin/bash

## Check flags
pash_output_time_flag=0
pash_execute_flag=1
for item in $@
do
    if [ "--output_time" == "$item" ]; then
        pash_output_time_flag=1
    fi

    if [ "--compile_optimize_only" == "$item" ] || [ "--compile_only" == "$item" ]; then
        pash_execute_flag=0
    fi
done

## The first argument contains the sequential script. Just running it should work for all tests.
>&2 echo "Sequential file in: $1 contains"
>&2 cat $1
## TODO: Don't write to /tmp_distr_output/0 anymore
source $1 > /tmp/distr_output/0

# pash_compiled_script_file=$(mktemp -u)
# python3.8 pash_runtime.py ${pash_compiled_script_file} $@

# ## Count the execution time and execute the compiled script
# pash_exec_time_start=$(date +"%s%N")
# if [ "$pash_execute_flag" -eq 1 ]; then
#     source ${pash_compiled_script_file}
# fi
# pash_exec_time_end=$(date +"%s%N")

# ## TODO: Maybe remove the temp file after execution

# ## We want the execution time in milliseconds
# if [ "$pash_output_time_flag" -eq 1 ]; then
#     pash_exec_time_ms=$(echo "scale = 3; ($pash_exec_time_end-$pash_exec_time_start)/1000000" | bc)
#     >&2 echo "Execution time: $pash_exec_time_ms  ms"
# fi
