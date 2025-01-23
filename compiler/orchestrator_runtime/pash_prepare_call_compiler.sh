#!/bin/bash

## OUTPUT: When it completes it sets "$pash_script_to_execute"

## Only needed for expansion
export pash_input_args=( "$@" )

## Move some pash env variables to local so that tests pass
tmp="$pash_disable_parallel_pipelines"
unset pash_disable_parallel_pipelines
pash_disable_parallel_pipelines="$tmp"

tmp="$pash_input_ir_file"
unset pash_input_ir_file
pash_input_ir_file="$tmp"

tmp="$pash_sequential_script_file"
unset pash_sequential_script_file
pash_sequential_script_file="$tmp"

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

if [[ "$daemon_response" == *"not all regions are parallelizable"* ]]; then
    pash_all_region_parallelizable=1 
else 
    pash_all_region_parallelizable=0 
fi

if [[ "$daemon_response" == *"OK:"* ]]; then
    pash_runtime_return_code=0
elif [ -z "$daemon_response" ]; then
    ## Trouble... Daemon crashed, rip
    pash_redir_output echo "$$: ERROR: (2) Daemon crashed!"
    exit 1
else
    pash_runtime_return_code=1
fi

# save IFS to restore after field splitting
[ -n "${IFS+x}" ] && saved_IFS=$IFS
unset IFS
# Get assigned process id
# We need to split the daemon response into elements of an array by
# shell's field splitting.
# shellcheck disable=SC2206
response_args=($daemon_response)
process_id=${response_args[1]}

[ -n "${saved_IFS+x}" ] && IFS="$saved_IFS"

pash_redir_output echo "$$: (2) Compiler exited with code: $pash_runtime_return_code"

## only when --assert_all_regions_parallellizable is used do we care about all regions being parallelizable
if [ "$pash_all_region_parallelizable" -ne 0 ] && [ "$pash_assert_all_regions_parallelizable_flag" -eq 1 ]; then
    pash_redir_output echo "$$: ERROR: (2) Compiler failed with error code because some regions were not parallelizable: $pash_all_region_parallelizable while assert_all_regions_parallelizable_flag was enabled! Exiting PaSh..."
    exit 1
fi

if [ "$pash_runtime_return_code" -ne 0 ] && [ "$pash_assert_all_regions_parallelizable_flag" -eq 1 ]; then
    pash_redir_output echo "$$: ERROR: (2) Compiler failed with error code: $pash_runtime_return_code while assert_all_regions_parallelizable_flag was enabled! Exiting PaSh..."
    exit 1
fi

## for pash_assert_compiler_success_flag, exit when return code is 0 (general exception caught) and not when all regions are parallelizable
if [ "$pash_runtime_return_code" -ne 0 ] && [ "$pash_all_region_parallelizable" -eq 0 ] && [ "$pash_assert_compiler_success_flag" -eq 1 ]; then
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

