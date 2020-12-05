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
pash_sequential_script_file=$1
# >&2 echo "Sequential file in: $1 contains"
# >&2 cat $1
# source $1

## Current plan: Optimistic sequential execution with eager in stdin stdout. 
## If the compiler succeeds in compiling (and improving) the script and the following two constraints hold:
##  (1) The outputs of the DFG are not appended
##  (2) All input and output files are simple files in the file system (no device files, fifos, etc)
## then we empty the output files, stop the original script, and execute the parallel script.
## (Note that we have to reroute what was in the stdin eager buffer to the parallel script.)
##
## Else if the original script finished execution before compilation succeeds, or if compilation fails:
## - We pipe stdout eager to stdout
## - We let everything finish.

## Assumptions/Constraints:
##
## For now we can assume that constraints (1) and (2) hold.
##
## TODO: A first TODO would be to check them in the compilation process
##
## TODO: An alternative TODO would be to let preprocessing give us information about them, allowing us to 
##       have a finer tuned execution plan depending on this information. For example, if we see that script
##       has append to some file we can be carefull and buffer its output using eager.

## Implementation proposal: Use a signal trap
##
## ISSUE: It seems that bash waits to handle a trap only after a command that it executes finished execution.
##        This means that we can't reliably use source and signal traps.
##
## POSSIBLE SOLUTION: 
##        At the moment we use source for two reasons.
##          1. To see changes in variables in the parent shell.
##          2. To see local variables in the sourced shell.
##        To address the above issue, we could avoid using source,
##        and actually run the original script in a subshell that first reads its parents
##        local variables and then exports back the variables.

## NOTE: The intuition about why quick-abort works is that if the compilation succeeds, then the
##       script is a DFG, meaning that we know exactly how it affects its environment after completing.
##       Therefore, we can go back and stop the already running script without risking unsafe behavior.

## TODO: We need to *at least* capture the stdout and stderr of the executed bash script 
##       if we ever want to have a change of running a compiled version in its place.
##       In order to do that, we can pipe its stdout and stderr to an eager command (so that it doesn't get lost)
##       and if we actually need it, we can pipe from the output of the eager command to stdout and stderr.
##
##       If it ends out that we don't need it we can just kill the command and dump its
##       outputs to /dev/null.

## TODO: If the output of the DFG is more than just stdout, then we actually need to first delete
##       what the outputs (since they will be rewritten from the DFG).

## ISSUE: What about if the DFG appends(!) to its output (instead of just writing). Then
##        we can't handle it correctly.

## ISSUE: We can't actually `source` the execution of the sequential if we later want to kill it.
##
## POSSIBLE SOLUTION: 
##        Maybe we can develop a signal handler in the PaSh runtime that receives a signal and starts executing the parallel version.
##        Then we can fork of a process that runs the python `pash_runtime` and if it is completed, it can send the signal to 
##        current shell.
## 
## TODO: If we choose to follow this solution we have to make sure that there is no race condition.

## TODO: We also want to avoid executing the compiled script if it doesn't contain any improvement.

## Prepare a file with all shell variables
pash_runtime_shell_variables_file=$(mktemp -u)
declare -p > $pash_runtime_shell_variables_file

pash_compiled_script_file=$(mktemp -u)
python3.8 pash_runtime.py ${pash_compiled_script_file} --var_file "${pash_runtime_shell_variables_file}" "${@:2}"
pash_runtime_return_code=$?

## Count the execution time and execute the compiled script
pash_exec_time_start=$(date +"%s%N")
if [ "$pash_execute_flag" -eq 1 ]; then
    ## If the compiler failed, we have to run the sequential
    if [ "$pash_runtime_return_code" -ne 0 ]; then
        source ${pash_sequential_script_file}
    else
        source ${pash_compiled_script_file}
    fi
fi
pash_exec_time_end=$(date +"%s%N")

## TODO: Maybe remove the temp file after execution

## We want the execution time in milliseconds
if [ "$pash_output_time_flag" -eq 1 ]; then
    pash_exec_time_ms=$(echo "scale = 3; ($pash_exec_time_end-$pash_exec_time_start)/1000000" | bc)
    >&2 echo "Execution time: $pash_exec_time_ms  ms"
fi
