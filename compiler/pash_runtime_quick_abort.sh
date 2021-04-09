#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

still_alive()
{
    jobs -p | tr '\n' ' '
}

log()
{
    pash_redir_output echo "$$: (QAbort) " "$@"
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

## TODO: Maybe the reroute needs to be put around (C) and not (D)

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
    seq_input_eager_pid=$!
    log "Spawned sequential input eager: $seq_input_eager_pid with:"
    log "  -- IN: $pash_tee_stdout1"
    log "  -- OUT: $pash_seq_eager_output"
    log "  -- INTERM: $pash_seq_eager_file"

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
    log "Sequential pid: $pash_seq_pid"

    ## (E) The sequential output eager
    pash_seq_eager2_output="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    mkfifo "$pash_seq_eager2_output"
    pash_seq_eager2_file="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    "$RUNTIME_DIR/../runtime/eager" "$pash_seq_output" "$pash_seq_eager2_output" "$pash_seq_eager2_file" &

    ## (F) Second eager
    pash_par_eager_output="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    mkfifo "$pash_par_eager_output"
    pash_par_eager_file="$($RUNTIME_DIR/pash_ptempfile_name.sh)"
    "$RUNTIME_DIR/../runtime/eager" "$pash_tee_stdout2" "$pash_par_eager_output" "$pash_par_eager_file" &
    par_eager_pid=$!
    log "Spawned parallel eager: $par_eager_pid with:"
    log "  -- IN: $pash_tee_stdout2"
    log "  -- OUT: $pash_par_eager_output"
    log "  -- INTERM: $pash_par_eager_file"

    ## Run the compiler
    pash_redir_all_output python3 "$RUNTIME_DIR/pash_runtime.py" ${pash_compiled_script_file} --var_file "${pash_runtime_shell_variables_file}" "${@:2}" &
    pash_compiler_pid=$!
    log "Compiler pid: $pash_compiler_pid"

    ## Wait until one of the two (original script, or compiler) die
    alive_pids=$(still_alive)
    log "Still alive: $alive_pids"
    while `list_include_item "$alive_pids" "$pash_seq_pid"` && `list_include_item "$alive_pids" "$pash_compiler_pid"` ; do
        ## Wait for either of the two to complete
        wait -n "$pash_seq_pid" "$pash_compiler_pid"
        completed_pid_status=$?
        log "Process exited with return code: $completed_pid_status"
        alive_pids=$(still_alive)
        log "Still alive: $alive_pids"
    done

    ## If the sequential is still alive we want to see if the compiler succeeded
    if `list_include_item "$alive_pids" "$pash_seq_pid"` ; then
        pash_runtime_return_code=$completed_pid_status
        log "Compilation was done first with return code: $pash_runtime_return_code"

        ## We only want to run the parallel if the compiler succeeded.
        ## TODO: Enable that
        if [ "$pash_runtime_return_code" -eq 0 ]; then

            ## TODO: We really need to kill the sequential (so that it stops writing to other outputs).
            ##       Actually we need to call it with reroute to dump its stdin to /dev/null and kill it.
            log "Killing sequential pid: $pash_seq_pid..."
            kill -n 9 "$pash_seq_pid" 2> /dev/null
            kill_status=$?
            wait "$pash_seq_pid" 2> /dev/null
            seq_exit_status=$?
            log "Sequential pid: $pash_seq_pid was killed successfully returning status $seq_exit_status."
            log "Still alive: $(still_alive)"

            ## If kill failed it means it was already completed, 
            ## and therefore we do not need to run the parallel.
            if true || [ "$kill_status" -eq 0 ]; then
                ## (2) Run the parallel
                log "Run parallel:"
                log "  -- Runtime vars: $pash_runtime_shell_variables_file"
                log "  -- Output vars: $pash_output_variables_file"
                log "  -- Output set: ${pash_output_set_file}"
                log "  -- Compiled script: ${pash_compiled_script_file}"
                log "  -- Input: $pash_par_eager_output"
                "$RUNTIME_DIR/pash_wrap_vars.sh" \
                    $pash_runtime_shell_variables_file \
                    $pash_output_variables_file \
                    ${pash_output_set_file} \
                    ${pash_compiled_script_file} \
                    < "$pash_par_eager_output"
                pash_runtime_final_status=$?
            fi
        else
            ## If the compiler failed we just wait until the sequential is done.

            ## (1) Redirect the seq output to stdout
            cat "$pash_seq_eager2_output" &
            log "STDOUT cat pid: $!"
            log "Still alive: $(still_alive)"

            log "Waiting for sequential: $pash_seq_pid"
            ## TODO: Do we need -n?
            # wait -n "$pash_seq_pid"
            wait "$pash_seq_pid"
            pash_runtime_final_status=$?
        fi
    else
        pash_runtime_final_status=$completed_pid_status
        log "Sequential was done first with return code: $pash_runtime_final_status"

        ## (1) Redirect the seq output to stdout
        cat "$pash_seq_eager2_output" &
        final_cat_pid=$!
        log "STDOUT cat pid: $final_cat_pid"

        ## If this fails (meaning that compilation is done) we do not care
        ## TODO: Do we actually need to kill compiler
        kill -9 "$pash_compiler_pid" 2> /dev/null
        wait "$pash_compiler_pid"  2> /dev/null
        log "Still alive: $(still_alive)"

        wait "$final_cat_pid"
    fi
    
    ## Return the exit code
    (exit "$pash_runtime_final_status")
fi
