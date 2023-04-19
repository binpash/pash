#!/bin/bash

## When it completes it sets "$pash_script_to_execute"

## Only needed for expansion
export pash_input_args=( "$@" )

## Save the shell variables to a file (necessary for expansion)
export pash_runtime_shell_variables_file="${PASH_TMP_PREFIX}/variables_$RANDOM$RANDOM$RANDOM"
source "$RUNTIME_DIR/pash_declare_vars.sh" "$pash_runtime_shell_variables_file"
pash_redir_output echo "$$: (1) Bash variables saved in: $pash_runtime_shell_variables_file"

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

## If the compiler failed or if we dry_run the compiler, we have to run the sequential
if [ "$pash_runtime_return_code" -ne 0 ] || [ "$pash_dry_run_compiler_flag" -eq 1 ]; then
    export pash_script_to_execute="${pash_sequential_script_file}"
else
    export pash_script_to_execute="${pash_compiled_script_file}"
fi

## Let daemon know that this region is done
function inform_daemon_exit () {
    ## Send to daemon
    msg="Exit:${process_id}"
    daemon_response=$(pash_communicate_daemon_just_send "$msg")
} 

