#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

still_alive()
{
    jobs -p | tr '\n' ' '
}

log_qa()
{
    log "(QAbort) " "$@"
}

log_time_from_qa()
{
    local start_time="$1"
    local prefix="$2"
    log_time_from "$start_time" "(QAbort) $prefix"
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

## This spawns a buffer command to buffer inputs and outputs
##
## It writes the pid to stdout
spawn_eager()
{
    local name=$1
    local input=$2
    local output=$3
    # local eager_file="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    ## Note: Using eager actually leads to some deadlock issues. It must have to do with eagers behavior when
    ##       its input or output closes.
    # "$RUNTIME_DIR/../runtime/eager" "$input" "$output" "$eager_file" </dev/null 1>/dev/null 2>/dev/null &
    "$RUNTIME_DIR/../runtime/dgsh_tee.sh" "$input" "$output" -I -f </dev/null 1>/dev/null 2>/dev/null &
    local eager_pid=$!
    log_qa "Spawned $name eager: $eager_pid with:"
    log_qa "  -- IN: $input"
    log_qa "  -- OUT: $output"
    # log_qa "  -- INTERM: $eager_file"
    echo "$eager_pid"
}

spawn_tee()
{
    local input=$1
    local output1=$2
    local output2=$3
    "$RUNTIME_DIR/../runtime/dgsh_tee" -i "$input" -o "$output1" -o "$output2" -f </dev/null 1>/dev/null 2>/dev/null &
    local tee_pid=$!
    log_qa "Spawned dgsh-tee: $tee_pid with:"
    log_qa "  -- IN: $input"
    log_qa "  -- OUT1: $output1"
    log_qa "  -- OUT1: $output2"
    echo "$tee_pid"
}

## Kills the process group that belongs to the given pgid
kill_pg()
{
    local pg_lead_pid=$1
    /bin/kill -15 "-${pg_lead_pid}" 2> /dev/null
}

## TODO: Make sure that this waits for all processes in the process group to finish executing.
wait_pg()
{
    local pg_lead_pid=$1
    wait "$pg_lead_pid" 2> /dev/null    
}

kill_wait_pg()
{
    kill_pg "$1"
    wait_pg "$1"
}

## Solution Schematic:
##
##  (A)      (B)      (C)       (D)        (E)
## stdin --- tee --- eager --- seq.sh --- eager --- OUT_SEQ
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

## TODO: Maybe the reroute needs to be put around (C) and not (D)

## TODO: Improve the happy path (very fast sequential) execution time 

## TODO: Use reroute around dgsh_tees to make sure that they do not use storage unnecessarily 
##       (if their later command is done).

## TODO: Catch SIGCHILD instead of waiting for everything to be set up before running.
##       This would require that we encapsulate all the sequential script execution in a different
##          script, to make sure that trapping SIGCHLD only happens if ony of the two big scripts exits.
##       Then we would have to have a big script for the compiler.
##       Only these two childs should be spawned, and then we can wait until any of them finishes.
##
##       Alternatively, we could implement this here in C, having finer-grained control over the concurrent processes.

if [ "$pash_execute_flag" -eq 1 ]; then
    ## We start with the sequential command so that it can finish first if it has no
    ## requirement for reading from stdin.
    pash_quick_abort_time_start=$(date +"%s%N")

    ## Setup all tempfiles at once
    pash_seq_eager_output="$(tempfile)"
    pash_seq_output="$(tempfile)"
    pash_seq_eager2_output="$(tempfile)"
    pash_tee_stdin="$(tempfile)"
    pash_tee_stdout1="$(tempfile)"
    pash_tee_stdout2="$(tempfile)"
    pash_par_eager_output="$(tempfile)"
    mkfifo "$pash_seq_eager_output" \
           "$pash_seq_output" \
           "$pash_seq_eager2_output" \
           "$pash_tee_stdin" \
           "$pash_tee_stdout1" \
           "$pash_tee_stdout2" \
           "$pash_par_eager_output"
    log_time_from_qa "$pash_quick_abort_time_start" "Tempfile preparation"

    pash_quick_abort_time_after_tempfiles=$(date +"%s%N")

    ## (D) Sequential command
    setsid "$RUNTIME_DIR/pash_wrap_vars.sh" \
        $pash_runtime_shell_variables_file \
        $pash_output_variables_file \
        ${pash_output_set_file} \
        ${pash_sequential_script_file} \
        > "$pash_seq_output" < "$pash_seq_eager_output" &
    pash_seq_pid=$!
    log_qa "Sequential pid: $pash_seq_pid"

    ## (E) The sequential output eager
    seq_output_eager_pid=$(spawn_eager "sequential output" "$pash_seq_output" "$pash_seq_eager2_output")

    ## (A) Redirect stdin to `tee`
    ## The redirections below are necessary to ensure that the background `cat` reads from stdin.
    { setsid cat > "$pash_tee_stdin" <&3 3<&- & } 3<&0
    pash_input_cat_pid=$!
    log_qa "Spawned input cat with pid: $pash_input_cat_pid"

    ## (B) A `tee` that duplicates input to both the sequential and parallel
    # tee_pid=$(spawn_tee "$pash_tee_stdin" "$pash_tee_stdout1" "$pash_tee_stdout2")
    tee "$pash_tee_stdout1" > "$pash_tee_stdout2" < "$pash_tee_stdin" &

    ## TODO: We probably don't need this.
    # (C) The sequential input eager
    seq_input_eager_pid=$(spawn_eager "sequential input" "$pash_tee_stdout1" "$pash_seq_eager_output")

    ## (F) Second eager
    par_eager_pid=$(spawn_eager "parallel input" "$pash_tee_stdout2" "$pash_par_eager_output")

    ## Run the compiler
    setsid python3 "$RUNTIME_DIR/pash_runtime.py" ${pash_compiled_script_file} --var_file "${pash_runtime_shell_variables_file}" "${@:2}" &
    pash_compiler_pid=$!
    log_qa "Compiler pid: $pash_compiler_pid"
    log_time_from_qa "$pash_quick_abort_time_after_tempfiles" "Phases (A-F)"

    ## Wait until one of the two (original script, or compiler) die
    pash_quick_abort_time_before_wait=$(date +"%s%N")
    ## TODO: Delete this loop if not necessary. It seems to be however...
    alive_pids=$(still_alive)
    log_qa "Still alive: $alive_pids"
    while `list_include_item "$alive_pids" "$pash_seq_pid"` && `list_include_item "$alive_pids" "$pash_compiler_pid"` ; do
        ## Wait for either of the two to complete
        wait -n "$pash_seq_pid" "$pash_compiler_pid"
        completed_pid_status=$?
        log_qa "Process exited with return code: $completed_pid_status"
        alive_pids=$(still_alive)
        log_qa "Still alive: $alive_pids"
    done
    log_time_from_qa "$pash_quick_abort_time_before_wait" "Waiting"

    ## If the sequential is still alive we want to see if the compiler succeeded
    if `list_include_item "$alive_pids" "$pash_seq_pid"` ; then
        pash_runtime_return_code=$completed_pid_status
        log_qa "Compilation was done first with return code: $pash_runtime_return_code"

        ## We only want to run the parallel if the compiler succeeded.
        if [ "$pash_runtime_return_code" -eq 0 ]; then

            ## TODO: Is this necessary
            ## Redirect the sequential output to /dev/null
            cat "$pash_seq_eager2_output" > /dev/null &
            seq_cat_pid=$!
            log_qa "seq to /dev/null cat pid: $seq_cat_pid"

            ## Kill the sequential process tree
            log_qa "Killing sequential pid: $pash_seq_pid..."
            kill_pg "$pash_seq_pid"
            kill_status=$?
            wait_pg "$pash_seq_pid"
            seq_exit_status=$?
            log_qa "Sequential pid: $pash_seq_pid was killed successfully returning status $seq_exit_status."

            ## If kill failed it means it was already completed, 
            ## and therefore we do not need to run the parallel.
            ##
            ## TOOD: Enable this optimization
            if true || [ "$kill_status" -eq 0 ]; then
                ## (2) Run the parallel
                log_qa "Run parallel:"
                log_qa "  -- Runtime vars: $pash_runtime_shell_variables_file"
                log_qa "  -- Output vars: $pash_output_variables_file"
                log_qa "  -- Output set: ${pash_output_set_file}"
                log_qa "  -- Compiled script: ${pash_compiled_script_file}"
                log_qa "  -- Input: $pash_par_eager_output"

                "$RUNTIME_DIR/pash_wrap_vars.sh" \
                    $pash_runtime_shell_variables_file \
                    $pash_output_variables_file \
                    ${pash_output_set_file} \
                    ${pash_compiled_script_file} \
                    < "$pash_par_eager_output" &
                ## Note: For some reason the above redirection used to create some issues,
                ##       but no more after we started using dgsh-tee

                pash_par_pid=$!
                log_qa "Parallel is running with pid: $pash_par_pid..."
                # strace -p $pash_par_pid 2>> $PASH_REDIR
                wait "$pash_par_pid"
                pash_runtime_final_status=$?
                log_qa "Parallel is done with status: $pash_runtime_final_status"
            else
                ## TODO: Handle that case properly by enabling the optimization above.
                log_qa "ERROR: Shouldn't have reached that"
                exit 1
            fi
        else
            ## If the compiler failed we just wait until the sequential is done.

            ## (1) Redirect the seq output to stdout
            cat "$pash_seq_eager2_output" &
            seq_output_cat_pid=$!
            log_qa "STDOUT cat pid: $seq_output_cat_pid"

            log_qa "Waiting for sequential: $pash_seq_pid"
            wait "$pash_seq_pid"
            pash_runtime_final_status=$?
            log_qa "DONE Sequential: $pash_seq_pid exited with status: $pash_runtime_final_status"

            ## TODO: It is not clear if we also need to wait for the output cat to end.
            log_qa "Waiting for sequential output cat: $seq_output_cat_pid"
            wait "$seq_output_cat_pid"
            log_qa "DONE Waiting for sequential output cat: $seq_output_cat_pid"

        fi
    else
        pash_runtime_final_status=$completed_pid_status
        log_qa "Sequential was done first with return code: $pash_runtime_final_status"

        ## (1) Redirect the seq output to stdout
        cat "$pash_seq_eager2_output" &
        final_cat_pid=$!
        log_qa "STDOUT cat pid: $final_cat_pid"

        ## We need to kill the compiler to not get delayed log output
        ## If this fails (meaning that compilation is done) we do not care
        kill_wait_pg "$pash_compiler_pid"

        wait "$final_cat_pid"
    fi

    pash_quick_abort_time_before_cleanup=$(date +"%s%N")
    ## TODO: Not clear if this is needed or if it doesn indeed kill all the
    ##       processes and cleans up everything properly
    ## Kill the input process
    log_qa "Killing the input cat process: $pash_input_cat_pid"
    kill_wait_pg "$pash_input_cat_pid"
    # kill -9 $pash_input_cat_pid 2> /dev/null
    # wait $pash_input_cat_pid 2> /dev/null
    log_qa "The input cat: $pash_input_cat_pid died!"
    

    ## TODO: This (and the above) should not be needed actually, everything should be already done due to
    ##       sequential and parallel both having exited.
    ## Kill every spawned process
    still_alive_pids="$(still_alive)"
    log_qa "Killing all the still alive: $still_alive_pids"
    kill -15 "$still_alive_pids" 2> /dev/null
    wait $still_alive_pids 2> /dev/null
    log_qa "All the alive pids died: $still_alive_pids"
    log_time_from_qa "$pash_quick_abort_time_before_cleanup" "Cleanup"

    ## Return the exit code
    (exit "$pash_runtime_final_status")
fi
