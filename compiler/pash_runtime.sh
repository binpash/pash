#!/bin/bash

##
## High level design. 
##
## (1) The `pash_runtime` should behave as a wrapper, saving all the necessary state:
##     - previous exit code
##     - previous set status
##     - previous variables
##     and then reverting to PaSh internal state
##
## (2) Then it should perform pash-internal work.
##
## (3) Then it should make sure to revert the exit code and `set` state to the saved values.
##
## (4) Then it should execute the inside script (either original or parallel)
##
## (5) Then it save all necessary state and revert to pash-internal state. 
##     (At the moment this happens automatically because the script is ran in a subshell.)
##
## (6) Then it should do all left pash internal work.
##
## (7) Before exiting it should revert all exit state.
##
## Visually:
##
## -- bash -- | -- pash --
##    ...     |
##      \----(1)----\
##            |     ...
##            |     (2)
##            |     ...
##      /----(3)----/
##    ...     |
##    (4)     |
##    ...     |
##
## (The rest of the steps happen only in debug mode)
##    ...
##      \----(5)----\
##            |     ...
##            |     (6)
##            |     ...
##      /----(7)----/
##    ...     |

##
## Necessary for bash:
## - Last exit code $?
## - set state $-
##

##
## (1)
##

## TODO: Replace these exports completely (and only leave the debug prints)
export pash_previous_exit_status="$PREVIOUS_SHELL_EC"
# export pash_input_args=( "${PREVIOUS_SHELL_ARGS[@]}" )
export pash_previous_set_status="$PREVIOUS_SET_STATUS"

pash_redir_output echo "$$: (1) Previous exit status: $pash_previous_exit_status"
pash_redir_output echo "$$: (1) Previous set state: $pash_previous_set_status"
pash_redir_output echo "$$: (1) Set state reverted to PaSh-internal set state: $-"


## Save the shell variables to a file
export pash_runtime_shell_variables_file="${PASH_TMP_PREFIX}/variables_$RANDOM$RANDOM$RANDOM"
source "$RUNTIME_DIR/pash_declare_vars.sh" "$pash_runtime_shell_variables_file"
pash_redir_output echo "$$: (1) Bash variables saved in: $pash_runtime_shell_variables_file"

##
## (2)
##

if [ "$pash_speculative_flag" -eq 1 ]; then
    ## For speculation, we don't want to invoke the compiler (for now) in (2),
    ## we just want to ask the scheduler in (3) to let us know when the df_region
    ## has finished executing and what is its exit code.

    ## The first argument is just the command id
    export pash_speculative_command_id=$1

    source "$RUNTIME_DIR/speculative/speculative_runtime.sh" "${pash_speculative_command_id}"

    ## TODO:
    ## 2. Check the flag in pash.py and if it is set, do the speculative transformation.
    

    ## TODO: (Future) We also want to let the scheduler know of any variable changes
    ## TODO: (Future) Check how we could support the steps (5), (6) with speculative and how to refactor this code the best way possible.
    ## TODO: (Future) We might not need all the set state and other config done in (1) and (3) for speculative
else

    ## The first argument contains the sequential script. Just running it should work for all tests.
    pash_sequential_script_file=$1

    ## The second argument SHOULD be the file that contains the IR to be compiled 
    pash_input_ir_file=$2

    ## The parallel script will be saved in the following file if compilation is successful.
    pash_compiled_script_file="${PASH_TMP_PREFIX}/pash_$RANDOM$RANDOM$RANDOM"


    ## TODO: Have a more proper communication protocol
    ## TODO: Make a proper client for the daemon
    pash_redir_output echo "$$: (2) Before asking the daemon for compilation..."
    ## Send and receive from daemon
    msg="Compile:${pash_compiled_script_file}| Variable File:${pash_runtime_shell_variables_file}| Input IR File:${pash_input_ir_file}"
    daemon_response=$(pash_communicate_daemon "$msg") # Blocking step, daemon will not send response until it's safe to continue

    if [[ "$daemon_response" == *"OK:"* ]]; then
        pash_runtime_return_code=0
    elif [ -z "$daemon_response" ]; then
        ## Trouble... Daemon crashed, rip
        pash_redir_output echo "$$: ERROR: (2) Daemon crashed!"
        exit 1
    else
        pash_runtime_return_code=1
    fi

    # Get assigned process id
    # We need to split the daemon response into elements of an array by
    # shell's field splitting.
    # shellcheck disable=SC2206
    response_args=($daemon_response)
    process_id=${response_args[1]}

    pash_redir_output echo "$$: (2) Compiler exited with code: $pash_runtime_return_code"
    if [ "$pash_runtime_return_code" -ne 0 ] && [ "$pash_assert_compiler_success_flag" -eq 1 ]; then
        pash_redir_output echo "$$: ERROR: (2) Compiler failed with error code: $pash_runtime_return_code while assert_compiler_success was enabled! Exiting PaSh..."
        exit 1
    fi

    # store functions for distributed execution
    if [ "$distributed_exec" -eq 1 ]; then
        declared_functions="${PASH_TMP_PREFIX}/pash_$RANDOM$RANDOM$RANDOM"
        declare -f > "$declared_functions"
        export declared_functions
    fi

    ##
    ## (3)
    ##

    ## Count the execution time
    pash_exec_time_start=$(date +"%s%N")

    ## If the compiler failed or if we dry_run the compiler, we have to run the sequential
    if [ "$pash_runtime_return_code" -ne 0 ] || [ "$pash_dry_run_compiler_flag" -eq 1 ]; then
        pash_script_to_execute="${pash_sequential_script_file}"
    else
        pash_script_to_execute="${pash_compiled_script_file}"
    fi

    ##
    ## (4)
    ##

    ## Let daemon know that this region is done
    function inform_daemon_exit () {
        ## Send to daemon
        msg="Exit:${process_id}"
        daemon_response=$(pash_communicate_daemon_just_send "$msg")
    } 

    function run_parallel() {
        trap inform_daemon_exit SIGTERM SIGINT EXIT
        source "$RUNTIME_DIR/pash_wrap_vars.sh" "$pash_script_to_execute"
        # internal_exec_status=$?
        ## TODO: Why do these need to happen here?
        # final_steps
        inform_daemon_exit
        ## TODO: Why do we need this to happen here?
        # (exit $internal_exec_status)
    }

    ## TODO: Delete this if useless
    ## We only want to execute (5) and (6) if we are in debug mode and it is not explicitly avoided
    # function final_steps() {
    #     if [ "$PASH_DEBUG_LEVEL" -ne 0 ] && [ "$pash_avoid_pash_runtime_completion_flag" -ne 1 ]; then
    #         ##
    #         ## (5)
    #         ##

    #         ## Prepare a file for the output shell variables to be saved in
    #         pash_output_var_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")
    #         # pash_redir_output echo "$$: Output vars: $pash_output_var_file"

    #         ## Prepare a file for the `set` state of the inner shell to be output
    #         pash_output_set_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")

    #         ## TODO: This can be turned to save_shell_state
    #         source "$RUNTIME_DIR/pash_runtime_shell_to_pash.sh" "$pash_output_var_file" "$pash_output_set_file"

    #         ##
    #         ## (6)
    #         ##
    #         source "$RUNTIME_DIR/pash_runtime_complete_execution.sh"
    #     fi
    # }

    ## TODO: Add a check that `set -e` is not on

    ## Check if there are traps set, and if so do not execute in parallel
    ##
    ## TODO: This might be an overkill but is conservative
    traps_set=$(trap)
    pash_redir_output echo "$$: (2) Traps set: $traps_set"
    # Don't fork if compilation failed. The script might have effects on the shell state.
    if [ "$pash_runtime_return_code" -ne 0 ] ||
        ## If parallel pipelines is not enabled we shouldn't fork 
        [ "$pash_parallel_pipelines" -eq 0 ] ||
        ## If parallel pipelines is explicitly disabled (e.g., due to context), no forking
        [ "$pash_disable_parallel_pipelines" -eq 1 ] ||
        ## If traps are set, no forking
        [ ! -z "$traps_set" ]; then
        # Early clean up in case the script effects shell like "break" or "exec"
        # This is safe because the script is run sequentially and the shell 
        # won't be able to move forward until this is finished

        ## Inform the daemon before we run
        ## TODO: Why not inform the daemon after? It should probably be after if we want to be correct
        inform_daemon_exit 
        source "$RUNTIME_DIR/pash_wrap_vars.sh" "$pash_script_to_execute"
        source "$RUNTIME_DIR/save_shell_state.sh"
        pash_runtime_final_status="$PREVIOUS_SHELL_EC"
        export pash_previous_set_status="$PREVIOUS_SET_STATUS"

        ## TODO: Move all this if/then/else in a different script
        if [ "$PASH_DEBUG_LEVEL" -ne 0 ] && [ "$pash_avoid_pash_runtime_completion_flag" -ne 1 ]; then
            ##
            ## (5)
            ##
            pash_redir_output echo "$$: (5) BaSh script exited with ec: $pash_runtime_final_status"
            pash_redir_output echo "$$: (5) Current BaSh shell: $pash_previous_set_status"
            pash_redir_output echo "$$: (5) Reverted to PaSh set state to: $-"

            ## Prepare a file for the output shell variables to be saved in
            pash_output_var_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")
            # pash_redir_output echo "$$: Output vars: $pash_output_var_file"
            source "$RUNTIME_DIR/pash_declare_vars.sh" "$pash_output_var_file"

            ## Prepare a file for the `set` state of the inner shell to be output
            pash_output_set_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")
            pash_redir_output echo "$$: (5) Writing current BaSh set state to: $pash_output_set_file"
            echo "$pash_previous_set_status" > "$pash_output_set_file"

            ##
            ## (6)
            ##
            source "$RUNTIME_DIR/pash_runtime_debug_complete_execution.sh"

            ## Restore the set state from a file because it has been rewritten by sourcing variables
            export pash_previous_set_status="$(cat "$pash_output_set_file")"
        fi
    else 
        # Should we redirect errors aswell?
        # TODO: capturing the return state here isn't completely correct. 
        run_parallel <&0 &
        ## Setting this to 0 since we can't capture this exit value
        pash_runtime_final_status=0
        pash_redir_output echo "$$: (2) Running pipeline..."

        ## The only thing we can recover here is the set state:
        ##   arguments, variables, and exit code cannot be returned
    fi
    ## Set the shell state before exiting
    pash_redir_output echo "$$: (7) Current PaSh set state: $-"
    source "$RUNTIME_DIR/pash_set_from_to.sh" "$-" "$pash_previous_set_status"
    pash_redir_output echo "$$: (7) Reverted to BaSh set state before exiting: $-"
fi
