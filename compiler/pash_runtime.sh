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
##     TODO: Figure out what could be different before (1), during (4), and after (7) 
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

## TODO: Make a list/properly define what needs to be saved at (1), (3), (5), (7)
##
## Necessary for pash:
## - PATH important for PaSh but might be changed in bash
## - IFS has to be kept default for PaSh to work
##
## Necessary for bash:
## - Last PID $! (TODO)
## - Last exit code $?
## - set state $-
## - File descriptors (TODO)
## - Loop state (?) Maybe `source` is adequate for this (TODO)
## - Traos (TODO)
##
## (maybe) TODO: After that, maybe we can create cleaner functions for (1), (3), (5), (7). 
##               E.g. we can have a correspondence between variable names and revert them using them 

##
## (1)
##

## Store the previous exit status to propagate to the compiler
## export pash_previous_exit_status=$?
## The assignment now happens outside
export pash_previous_exit_status

## Store the current `set` status to pash to the inside script 
export pash_previous_set_status=$-

pash_redir_output echo "$$: (1) Previous exit status: $pash_previous_exit_status"
pash_redir_output echo "$$: (1) Previous set state: $pash_previous_set_status"

## Prepare a file with all shell variables
##
## This is only needed by PaSh to expand.
##
## TODO: Maybe we can get rid of it since PaSh has access to the environment anyway?
## TODO: Remove this call to pash_ptempfile_name.sh. Actually remove this file in general.
##       PaSh should only generate temp files using $RANDOM$RANDOM$RANDOM
# pash_runtime_shell_variables_file="$($RUNTIME_DIR/pash_ptempfile_name.sh $distro)"
pash_runtime_shell_variables_file="${PASH_TMP_PREFIX}/pash_$RANDOM$RANDOM$RANDOM"
source "$RUNTIME_DIR/pash_declare_vars.sh" "$pash_runtime_shell_variables_file"
pash_redir_output echo "$$: (1) Bash variables saved in: $pash_runtime_shell_variables_file"

## Abort script if variable is unset
pash_default_set_state="huB"

## Revert the `set` state to not have spurious failures 
pash_redir_output echo "$$: (1) Bash set state at start of execution: $pash_previous_set_status"
source "$RUNTIME_DIR/pash_set_from_to.sh" "$pash_previous_set_status" "$pash_default_set_state"
pash_redir_output echo "$$: (1) Set state reverted to PaSh-internal set state: $-"

##
## (2)
##

## The first argument contains the sequential script. Just running it should work for all tests.
pash_sequential_script_file=$1

## The second argument SHOULD be the file that contains the IR to be compiled 
pash_input_ir_file=$2

## The parallel script will be saved in the following file if compilation is successful.
# pash_compiled_script_file="$($RUNTIME_DIR/pash_ptempfile_name.sh $distro)"
pash_compiled_script_file="${PASH_TMP_PREFIX}/pash_$RANDOM$RANDOM$RANDOM"


if [ "$pash_speculation_flag" -eq 1 ]; then
    ## Count the execution time
    pash_exec_time_start=$(date +"%s%N")
    source "$RUNTIME_DIR/pash_runtime_quick_abort.sh"
    pash_runtime_final_status=$?
    ## For now this will fail!!!
    exit 1
else

    if [ "$pash_daemon" -eq 1 ]; then
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
    else
        pash_redir_all_output_always_execute python3 -S "$RUNTIME_DIR//pash_runtime.py" --var_file "${pash_runtime_shell_variables_file}" "${pash_compiled_script_file}" "${pash_input_ir_file}" "$@"
        pash_runtime_return_code=$?
    fi

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
    
    # ##
    # ## (4)
    # ##

    ## TODO: It might make sense to move these functions in pash_init_setup to avoid the cost of redefining them here.
    function clean_up () {
        if [ "$pash_daemon" -eq 1 ]; then
            if [ "$parallel_script_time_start" == "None" ] || [ "$pash_profile_driven_flag" -eq 0 ]; then
                exec_time=""
            else
                parallel_script_time_end=$(date +"%s%N")
                parallel_script_time_ms=$(echo "scale = 3; ($parallel_script_time_end-$parallel_script_time_start)/1000000" | bc)
                pash_redir_output echo " --- --- Execution time: $parallel_script_time_ms  ms"
                exec_time=$parallel_script_time_ms
            fi
            ## Send to daemon
            msg="Exit:${process_id}|Time:$exec_time"
            daemon_response=$(pash_communicate_daemon_just_send "$msg")
        fi
    } 

    function run_parallel() {
        trap clean_up SIGTERM SIGINT EXIT
        if [ "$pash_profile_driven_flag" -eq 1 ]; then
            parallel_script_time_start=$(date +"%s%N")
        fi
        source "$RUNTIME_DIR/pash_wrap_vars.sh" "$pash_script_to_execute"
        internal_exec_status=$?
        final_steps
        clean_up
        (exit $internal_exec_status)
    }

    ## We only want to execute (5) and (6) if we are in debug mode and it is not explicitly avoided
    function final_steps() {
        if [ "$PASH_DEBUG_LEVEL" -ne 0 ] && [ "$pash_avoid_pash_runtime_completion_flag" -ne 1 ]; then
            ##
            ## (5)
            ##

            ## Prepare a file for the output shell variables to be saved in
            pash_output_var_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")
            # pash_redir_output echo "$$: Output vars: $pash_output_var_file"

            ## Prepare a file for the `set` state of the inner shell to be output
            pash_output_set_file=$("$RUNTIME_DIR/pash_ptempfile_name.sh" "$distro")

            source "$RUNTIME_DIR/pash_runtime_shell_to_pash.sh" "$pash_output_var_file" "$pash_output_set_file"

            ##
            ## (6)
            ##
            source "$RUNTIME_DIR/pash_runtime_complete_execution.sh"
        fi
    }

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
       [ ! -z "$traps_set" ] ||
       [ "$pash_daemon" -eq 0 ]; then
        # Early clean up in case the script effects shell like "break" or "exec"
        # This is safe because the script is run sequentially and the shell 
        # won't be able to move forward until this is finished

        ## Needed to clear up any past script time start execution times.        
        parallel_script_time_start=None
        clean_up 
        source "$RUNTIME_DIR/pash_wrap_vars.sh" "$pash_script_to_execute"
        pash_runtime_final_status=$?
        final_steps
    else 
        # Should we redirect errors aswell?
        # TODO: capturing the return state here isn't completely correct. 
        # Might need more complex design if this end up being a problem
        run_parallel <&0 &
        pash_runtime_final_status=$?
        pash_redir_output echo "$$: (2) Running pipeline"

        ## Here we need to also revert the state back to bash state 
        ## since run_parallel will do that in a separate shell
        ##
        ## This happens right before we exit from pash_runtime!

        ## Recover the `set` state of the previous shell
        # pash_redir_output echo "$$: (3) Previous BaSh set state: $pash_previous_set_status"
        # pash_redir_output echo "$$: (3) PaSh-internal set state of current shell: $-"
        pash_current_set_state=$-
        source "$RUNTIME_DIR/pash_set_from_to.sh" "$pash_current_set_state" "$pash_previous_set_status"
        pash_redir_output echo "$$: (5) Reverted to BaSh set state: $-"

        ## TODO: This might not be necessary
        ## Recover the input arguments of the previous script
        ## Note: We don't need to care about wrap_vars arguments because we have stored all of them already.
        #
        # This variable stores arguments as a space-separated stirng, so we
        # need to unquote it and to split it into multiple strings by shell's
        # field splitting.
        # shellcheck disable=SC2086
        set -- $pash_input_args
        pash_redir_output echo "$$: (5) Reverted to BaSh input arguments: $@"

        ## TODO: We probably need to exit with the exit code here or something!
    fi
fi

