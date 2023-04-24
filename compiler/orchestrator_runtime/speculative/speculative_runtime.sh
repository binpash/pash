#!/bin/bash


## TODO: Define the client in pash_spec_init_setup (which should be sourced by pash_init_setup)

export pash_speculative_command_id=$1

pash_redir_output echo "$$: (2) Before asking the scheduler for cmd: ${pash_speculative_command_id} exit code..."

## TODO: Correctly save variables
## Save the shell variables to a file (necessary for expansion)
export pash_runtime_shell_variables_file="${PASH_TMP_PREFIX}/variables_$RANDOM$RANDOM$RANDOM"
source "$RUNTIME_DIR/pash_declare_vars.sh" "$pash_runtime_shell_variables_file"
pash_redir_output echo "$$: (1) Bash variables saved in: $pash_runtime_shell_variables_file"


## Send and receive from daemon
msg="Wait:${pash_speculative_command_id}"
daemon_response=$(pash_spec_communicate_scheduler "$msg") # Blocking step, daemon will not send response until it's safe to continue

## Receive an exit code
if [[ "$daemon_response" == *"OK:"* ]]; then
    # shellcheck disable=SC2206
    response_args=($daemon_response)
    pash_redir_output echo "$$: (2) Scheduler responded: $daemon_response"
    cmd_exit_code=${response_args[1]}
    output_variable_file=${response_args[2]}
elif [ -z "$daemon_response" ]; then
    ## Trouble... Daemon crashed, rip
    pash_redir_output echo "$$: ERROR: (2) Scheduler crashed!"
    exit 1
else
    pash_redir_output echo "$$: ERROR: (2) Scheduler responded garbage ${daemon_response}!"
    exit 1
fi


pash_redir_output echo "$$: (2) Scheduler returned exit code: ${cmd_exit_code} for cmd with id: ${pash_speculative_command_id}."


pash_runtime_final_status=${cmd_exit_code}

## TODO: Restore the variables
pash_redir_output echo "$$: (2) Recovering script variables from: $output_variable_file"
source "$RUNTIME_DIR/pash_source_declare_vars.sh" "$output_variable_file"

## TODO: Also need to use wrap_vars maybe to `set` properly etc
