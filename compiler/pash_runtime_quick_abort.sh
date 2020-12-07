#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

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

## Execute both the sequential and the compiler in parallel and only run the parallel
## compiled script if the compiler finished before the sequential.

still_alive()
{
    jobs -p | tr '\n' ' '
}

# Taken from: https://stackoverflow.com/a/20473191
# list_include_item "10 11 12" "2"
function list_include_item {
  local list="$1"
  local item="$2"
  if [[ $list =~ (^|[[:space:]])"$item"($|[[:space:]]) ]] ; then
    # yes, list include item
    result=0
  else
    result=1
  fi
  return $result
}


if [ "$pash_execute_flag" -eq 1 ]; then
    ## Note: First solution does not contain eager in stdout, therefore not being interactive.
    ## TODO: (Optimization) Fix the issue with eager in the output to make it interactive

    ## TODO: Fix the issue with the input. We have to first redirect input through an eager
    ##       so that it is saved and then redirected to the parallel one.

    ## Before running the original script we need to redirect its input and output to 
    ## eager and pipes.
    pash_stdin_redir1=$(mktemp -u)
    pash_stdin_redir2=$(mktemp -u)
    # pash_stdout_redir1=$(mktemp -u)
    pash_stdout_redir2=$(mktemp -u)
    pash_stdin_eager_file=$(mktemp -u)
    >&2 echo "eager intermediate file: $pash_stdin_eager_file"
    # pash_stdout_eager_file=$(mktemp -u)
    mkfifo $pash_stdin_redir1 $pash_stdin_redir2
    # mkfifo $pash_stdin_redir1 $pash_stdin_redir2 $pash_stdout_redir1 $pash_stdout_redir2


    ## TODO: Find the eager directory correctly
    "$RUNTIME_DIR/../evaluation/tools/eager" "$pash_stdin_redir1" "$pash_stdin_redir2" "$pash_stdin_eager_file" &
    >&2 echo "STDIN eager pid: $!"
    # ../evaluation/tools/eager "$pash_stdout_redir1" "$pash_stdout_redir2" "$pash_stdout_eager_file" &
    # >&2 echo "STDOUT eager pid: $!"
    ## The redirections below are necessary to ensure that the background `cat` reads from stdin.
    { cat > "$pash_stdin_redir1" <&3 3<&- & } 3<&0
    >&2 echo "STDIN cat pid: $!"
    ## Note: We don't connect stdout_redir2 yet, since it has to be bufferred for correctness. 

    ## Run the original script
    # "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_sequential_script_file} > "$pash_stdout_redir1" < "$pash_stdin_redir2" &
    "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_sequential_script_file} > "$pash_stdout_redir2" < "$pash_stdin_redir2" &
    # "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_sequential_script_file} > "$pash_stdout_redir2" &
    pash_seq_pid=$!
    >&2 echo "Sequential pid: $pash_seq_pid"

    ## Run the compiler
    python3.8 pash_runtime.py ${pash_compiled_script_file} --var_file "${pash_runtime_shell_variables_file}" "${@:2}" &
    pash_compiler_pid=$!
    >&2 echo "Compiler pid: $pash_compiler_pid"

    
    ## Wait until one of the two (original script, or compiler) die
    alive_pids=$(still_alive)
    >&2 echo "Still alive: $alive_pids"
    while `list_include_item "$alive_pids" "$pash_seq_pid"` && `list_include_item "$alive_pids" "$pash_compiler_pid"` ; do
        ## Wait for either of the two to complete
        wait -n "$pash_seq_pid" "$pash_compiler_pid"
        completed_pid_status=$?
        alive_pids=$(still_alive)
        >&2 echo "Still alive: $alive_pids"
    done

    ## If the sequential is still alive we want to see if the compiler succeeded
    if `list_include_item "$alive_pids" "$pash_seq_pid"` ; then
    # if [ "$pash_seq_pid" -eq "$alive_pids" ]; then
        pash_runtime_return_code=$completed_pid_status
        >&2 echo "Compilation was done first with return code: $pash_runtime_return_code"

        ## We only want to run the parallel if the compiler succeeded.
        if [ "$pash_runtime_return_code" -eq 0 ]; then
            kill -n 9 "$pash_seq_pid" 2> /dev/null
            kill_status=$?
            wait "$pash_seq_pid" 2> /dev/null
            pash_runtime_final_status=$?
            >&2 echo "Still alive: $(still_alive)"

            ## If kill failed it means it was already completed, 
            ## and therefore we do not need to run the parallel.
            if [ "$kill_status" -eq 0 ]; then
                >&2 echo "Run parallel"
                ## TODO: Find the outputs/inputs of the DFG and make sure that the outputs 
                ##       are clean and normal files (and the inputs are normal files)

                ## TODO: Redirect stdin/stdout so that the parallel one gets them from the start
                ##       and so that there are no duplicate entries in the stdout.
                ## TODO: This seems like a non-trivial solution since it requires keeping stdin open, 
                ##       while "restarting" the eager, making it send its output to a new pipe, and
                ##       then first reading all from its intermediate file. 
                "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_compiled_script_file}
                pash_runtime_final_status=$?
            fi
        else
            ## If the compiler failed we just wait until the sequential is done.

            wait -n "$pash_seq_pid"

            ## TODO: Redirect eagers to stdin + stdout
            cat "$pash_stdout_redir2" &
            >&2 echo "STDOUT cat pid: $!"
            >&2 echo "Still alive: $(still_alive)"

            pash_runtime_final_status=$?
        fi
    else
        >&2 echo "Sequential was done first!"

        ## TODO: Redirect eagers to stdin + stdout
        cat "$pash_stdout_redir2" &
        >&2 echo "STDOUT cat pid: $!"

        ## If this fails (meaning that compilation is done) we do not care
        kill -n 9 "$pash_compiler_pid" 2> /dev/null
        wait -n "$pash_compiler_pid"  2> /dev/null
        pash_runtime_final_status=$completed_pid_status
        >&2 echo "Still alive: $(still_alive)"
    fi
fi