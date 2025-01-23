#!/bin/bash

## Assumes the following two variables are set:
## pash_sequential_script_file: the sequential script. Just running it should work for all tests.
## pash_input_ir_file: the file that contains the IR to be compiled 

## When called by spec, assumes this variable is set:
## pash_spec_command_id: the node id for the specific command


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
## (1)
##

## TODO: Make a shell script that is called wrap and
##       takes three arguments (pre_script, script, post_script)
##       and wraps a shell script file with a pre and a post code.
##
##       This wrap function should perform only the necessary shell setting,
##       exit code, and $!, transfer. Arguments should ideally be transfered too but
##       I don't know of a good way to do it for an external script from an internal one.
##
##       Maybe it can happen with eval
##
## The challenging aspect is how to make this work for the parallel pipelines

## First save the state of the shell
source "$RUNTIME_DIR/save_shell_state.sh"
## Rename variables to pash specific names
export pash_previous_exit_status="$PREVIOUS_SHELL_EC"
export pash_previous_set_status="$PREVIOUS_SET_STATUS"

pash_redir_output echo "$$: (1) Previous exit status: $pash_previous_exit_status"
pash_redir_output echo "$$: (1) Previous set state: $pash_previous_set_status"
pash_redir_output echo "$$: (1) Set state reverted to PaSh-internal set state: $-"

##
## (2)
##

if [ "$pash_speculative_flag" -eq 1 ]; then
    ## For speculation, we don't want to invoke the compiler (for now) in (2),
    ## we just want to ask the scheduler in (3) to let us know when the df_region
    ## has finished executing and what is its exit code.

    ## TODO: We probably need to make this a local variable so that POSIX tests pass
    export pash_speculative_command_id=$pash_spec_command_id

    source "$RUNTIME_DIR/speculative/speculative_runtime.sh"

    ## TODO:
    ## 2. Check the flag in pash.py and if it is set, do the speculative transformation.
    

    ## TODO: (Future) We also want to let the scheduler know of any variable changes
    ## TODO: (Future) Check how we could support the steps (5), (6) with speculative and how to refactor this code the best way possible.
    ## TODO: (Future) We might not need all the set state and other config done in (1) and (3) for speculative
else
    ## Invoke the compiler and make any necessary preparations
    source "$RUNTIME_DIR/pash_prepare_call_compiler.sh"

    ## Check if there are traps set, and if so do not execute in parallel
    ## TODO: This might be an overkill but is conservative
    traps_set=$(trap)
    pash_redir_output echo "$$: (2) Traps set: $traps_set"
    # Don't fork if compilation failed. The script might have effects on the shell state.
    if [ "$pash_runtime_return_code" -ne 0 ] ||
        ## If parallel pipelines is disabled using a flag we shouldn't fork 
        [ "$pash_no_parallel_pipelines" -eq 1 ] ||
        ## If parallel pipelines is explicitly disabled (e.g., due to context), no forking
        [ "$pash_disable_parallel_pipelines" -eq 1 ] ||
        ## If traps are set, no forking
        [ ! -z "$traps_set" ]; then
        # Early clean up in case the script effects shell like "break" or "exec"
        # This is safe because the script is run sequentially and the shell 
        # won't be able to move forward until this is finished

        ## Inform the daemon (this happens before because otherwise when set -e is set we don't send the inform exit)
        ## However, this doesn't allow the compiler to get the proper execution time for a command
        ## TODO: Properly set and restore traps and then move inform afterwards
        ##       First make a test that has set traps and set -e to exit (check set-e.sh)
        ##
        ## TODO: Also inform the daemon that the timing does not work now so that it
        ##       doesn't measure time for profile driven optimizations.
        inform_daemon_exit 
        # echo $traps_set

        ## Run the script
        export SCRIPT_TO_EXECUTE="$pash_script_to_execute"
        source "$RUNTIME_DIR/pash_restore_state_and_execute.sh"
        ## Save the state after execution
        source "$RUNTIME_DIR/save_shell_state.sh"
        ## We don't need to save the arguments because they are already set
        pash_runtime_final_status="$PREVIOUS_SHELL_EC"
        export pash_previous_set_status="$PREVIOUS_SET_STATUS"

        pash_redir_output echo "$$: (5) BaSh script exited with ec: $pash_runtime_final_status"
    else 
        function run_parallel() {
            trap inform_daemon_exit SIGTERM SIGINT EXIT
            export SCRIPT_TO_EXECUTE="$pash_script_to_execute"
            source "$RUNTIME_DIR/pash_restore_state_and_execute.sh"
            inform_daemon_exit
        }
        # Should we redirect errors aswell?
        # TODO: capturing the return state here isn't completely correct. 
        run_parallel "$@" <&0 &
        ## Setting this to 0 since we can't capture this exit value
        pash_runtime_final_status=0
        pash_redir_output echo "$$: (2) Running pipeline..."

        ## The only thing we can recover here is the set state:
        ##  - arguments and variables are not modified since it is run in parallel and thus is pure
        ##  - exit code cannot be returned
    fi
    ## Set the shell state before exiting
    pash_redir_output echo "$$: (7) Current PaSh set state: $-"
    source "$RUNTIME_DIR/pash_set_from_to.sh" "$-" "$pash_previous_set_status"
    pash_redir_output echo "$$: (7) Reverted to BaSh set state before exiting: $-"
    ## Set the exit code
    (exit "$pash_runtime_final_status")
fi
