#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

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


## Solution Schematic:
##
##  (A)      (B)      (C)       (D)                (E)
## stdin --- tee --- eager --- reroute seq.sh --- eager --- OUT_SEQ
##            \      (F)
##             \--- eager --- PAR_IN
##
## (1) If compiler fails, or sequential is done executing:
##     - cat OUT_SEQ > stdout
##
## (2) If compiler succeeds:
##     - USR1 to reroute so that it redirects to /dev/null
##     - PAR_IN redirect to par stdin.
##
## Simplifying assumptions:
## - Not worrying about stderr
## - Not worrying about other inputs at the moment (assuming they are files if compiler succeeds)
## - Not worrying about other outputs 
##   + assuming that the parallel implementation will overwrite them
##   + Assuming that the DFG outputs are not appended
##
## TODO: A first TODO would be to check them in the compilation process
##
## TODO: An alternative TODO would be to let preprocessing give us information about them, allowing us to 
##       have a finer tuned execution plan depending on this information. For example, if we see that script
##       has append to some file we can be carefull and buffer its output using eager.

## NOTE: The intuition about why quick-abort works is that if the compilation succeeds, then the
##       script is a DFG, meaning that we know exactly how it affects its environment after completing.
##       Therefore, we can go back and stop the already running script without risking unsafe behavior.

## TODO: We also want to avoid executing the compiled script if it doesn't contain any improvement.

if [ "$pash_execute_flag" -eq 1 ]; then
    # set -x
    ## (A) Redirect stdin to `tee`
    pash_tee_stdin="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    mkfifo "$pash_tee_stdin"
    ## The redirections below are necessary to ensure that the background `cat` reads from stdin.
    { cat > "$pash_tee_stdin" <&3 3<&- & } 3<&0

    ## (B) A `tee` that duplicates input to both the sequential and parallel
    pash_tee_stdout1="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    pash_tee_stdout2="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    mkfifo "$pash_tee_stdout1" "$pash_tee_stdout2"
    tee "$pash_tee_stdout1" > "$pash_tee_stdout2" < "$pash_tee_stdin" &

    ## (C) The sequential input eager
    pash_seq_eager_output="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    mkfifo "$pash_seq_eager_output"
    pash_seq_eager_file="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    "$RUNTIME_DIR/../runtime/eager" "$pash_tee_stdout1" "$pash_seq_eager_output" "$pash_seq_eager_file" &

    ## (D) Sequential command
    pash_seq_output="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    mkfifo "$pash_seq_output"
    "$RUNTIME_DIR/pash_wrap_vars.sh" \
        $pash_runtime_shell_variables_file \
        $pash_output_variables_file \
        ${pash_output_set_file} \
        ${pash_sequential_script_file} \
        > "$pash_seq_output" < "$pash_seq_eager_output" &
    pash_seq_pid=$!
    pash_redir_output echo "Sequential pid: $pash_seq_pid"

    ## (F) Second eager
    cat "$pash_tee_stdout2" > /dev/null &

    cat "$pash_seq_output"

    wait $pash_seq_pid
    pash_seq_status=$?
    (exit "$pash_seq_status")


    # ## Before running the original script we need to redirect its input and output to 
    # ## eager and pipes.
    # pash_stdin_redir1="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    # pash_stdin_redir2="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    # # pash_stdout_redir1="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    # pash_stdout_redir2="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    # pash_stdin_eager_file="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    # pash_redir_output echo "eager intermediate file: $pash_stdin_eager_file"
    # # pash_stdout_eager_file="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    # mkfifo $pash_stdin_redir1 $pash_stdin_redir2
    # # mkfifo $pash_stdin_redir1 $pash_stdin_redir2 $pash_stdout_redir1 $pash_stdout_redir2


    # ## TODO: Find the eager directory correctly
    # "$RUNTIME_DIR/../evaluation/tools/eager" "$pash_stdin_redir1" "$pash_stdin_redir2" "$pash_stdin_eager_file" &
    # pash_redir_output echo "STDIN eager pid: $!"
    # # ../evaluation/tools/eager "$pash_stdout_redir1" "$pash_stdout_redir2" "$pash_stdout_eager_file" &
    # # pash_redir_output echo "STDOUT eager pid: $!"
    
    # pash_redir_output echo "STDIN cat pid: $!"
    # ## Note: We don't connect stdout_redir2 yet, since it has to be bufferred for correctness. 

    # ## Run the original script
    # # "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_output_set_file} ${pash_sequential_script_file} > "$pash_stdout_redir1" < "$pash_stdin_redir2" &
    # "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_output_set_file} ${pash_sequential_script_file} > "$pash_stdout_redir2" < "$pash_stdin_redir2" &
    # # "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_output_set_file} ${pash_sequential_script_file} > "$pash_stdout_redir2" &
    # pash_seq_pid=$!
    # pash_redir_output echo "Sequential pid: $pash_seq_pid"

    # ## Run the compiler
    # pash_redir_all_output python3 pash_runtime.py ${pash_compiled_script_file} --var_file "${pash_runtime_shell_variables_file}" "${@:2}" &
    # pash_compiler_pid=$!
    # pash_redir_output echo "Compiler pid: $pash_compiler_pid"

    
    # ## Wait until one of the two (original script, or compiler) die
    # alive_pids=$(still_alive)
    # pash_redir_output echo "Still alive: $alive_pids"
    # while `list_include_item "$alive_pids" "$pash_seq_pid"` && `list_include_item "$alive_pids" "$pash_compiler_pid"` ; do
    #     ## Wait for either of the two to complete
    #     wait -n "$pash_seq_pid" "$pash_compiler_pid"
    #     completed_pid_status=$?
    #     alive_pids=$(still_alive)
    #     pash_redir_output echo "Still alive: $alive_pids"
    # done

    # ## If the sequential is still alive we want to see if the compiler succeeded
    # if `list_include_item "$alive_pids" "$pash_seq_pid"` ; then
    # # if [ "$pash_seq_pid" -eq "$alive_pids" ]; then
    #     pash_runtime_return_code=$completed_pid_status
    #     pash_redir_output echo "Compilation was done first with return code: $pash_runtime_return_code"

    #     ## We only want to run the parallel if the compiler succeeded.
    #     if [ "$pash_runtime_return_code" -eq 0 ]; then
    #         kill -n 9 "$pash_seq_pid" 2> /dev/null
    #         kill_status=$?
    #         wait "$pash_seq_pid" 2> /dev/null
    #         pash_runtime_final_status=$?
    #         pash_redir_output echo "Still alive: $(still_alive)"

    #         ## If kill failed it means it was already completed, 
    #         ## and therefore we do not need to run the parallel.
    #         if [ "$kill_status" -eq 0 ]; then
    #             pash_redir_output echo "Run parallel"
    #             ## TODO: Find the outputs/inputs of the DFG and make sure that the outputs 
    #             ##       are clean and normal files (and the inputs are normal files)

    #             ## TODO: Redirect stdin/stdout so that the parallel one gets them from the start
    #             ##       and so that there are no duplicate entries in the stdout.
    #             ## TODO: This seems like a non-trivial solution since it requires keeping stdin open, 
    #             ##       while "restarting" the eager, making it send its output to a new pipe, and
    #             ##       then first reading all from its intermediate file. 
    #             "$RUNTIME_DIR/pash_wrap_vars.sh" $pash_runtime_shell_variables_file $pash_output_variables_file ${pash_output_set_file} ${pash_compiled_script_file}
    #             pash_runtime_final_status=$?
    #         fi
    #     else
    #         ## If the compiler failed we just wait until the sequential is done.

    #         wait -n "$pash_seq_pid"

    #         ## TODO: Redirect eagers to stdin + stdout
    #         cat "$pash_stdout_redir2" &
    #         pash_redir_output echo "STDOUT cat pid: $!"
    #         pash_redir_output echo "Still alive: $(still_alive)"

    #         pash_runtime_final_status=$?
    #     fi
    # else
    #     pash_redir_output echo "Sequential was done first!"

    #     ## TODO: Redirect eagers to stdin + stdout
    #     cat "$pash_stdout_redir2" &
    #     pash_redir_output echo "STDOUT cat pid: $!"

    #     ## If this fails (meaning that compilation is done) we do not care
    #     kill -n 9 "$pash_compiler_pid" 2> /dev/null
    #     wait -n "$pash_compiler_pid"  2> /dev/null
    #     pash_runtime_final_status=$completed_pid_status
    #     pash_redir_output echo "Still alive: $(still_alive)"
    # fi
fi
