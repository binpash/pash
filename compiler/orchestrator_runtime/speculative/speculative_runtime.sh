#!/bin/bash


## TODO: Define the client in pash_spec_init_setup (which should be sourced by pash_init_setup)

pash_redir_output echo "$$: (2) Before asking the scheduler for cmd: ${pash_speculative_command_id} exit code..."

## TODO: Correctly save variables
## Save the shell variables to a file (necessary for expansion)
export pash_runtime_shell_variables_file="${PASH_TMP_PREFIX}/variables_$RANDOM$RANDOM$RANDOM"
unset cmd_exit_code
unset output_variable_file
unset stdout_file
set +u
if [ -z "${PASH_OLD_IFS+x}" ]; then
    unset IFS
else
    IFS="$PASH_OLD_IFS"
fi
source "$RUNTIME_DIR/pash_declare_vars.sh" "$pash_runtime_shell_variables_file"
IFS=$' \t\n'
pash_redir_output echo "$$: (1) Bash variables saved in: $pash_runtime_shell_variables_file"

## TODO: We want to send the environment to the scheduler.
##       Once the scheduler determines if there are environment changes, it can then
##       decide to rerun or not the speculated commands with the new environment.


## Determine all current loop iterations and send them to the scheduler
pash_loop_iter_counters=${pash_loop_iters:-None}
pash_redir_output echo "$$: Loop node iteration counters: $pash_loop_iter_counters"

## Send and receive from daemon
msg="Wait:${pash_speculative_command_id}|Loop iters:${pash_loop_iter_counters}|Variables file:${pash_runtime_shell_variables_file}"
daemon_response=$(pash_spec_communicate_scheduler "$msg") # Blocking step, daemon will not send response until it's safe to continue

## Receive an exit code
if [[ "$daemon_response" == *"OK:"* ]]; then
    # shellcheck disable=SC2206
    response_args=($daemon_response)
    pash_redir_output echo "$$: (2) Scheduler responded: $daemon_response"
    
    cmd_exit_code=${response_args[1]}
    output_variable_file=${response_args[2]}
    stdout_file=${response_args[3]}

    pash_redir_output echo "$$: (2) Recovering stdout from: $stdout_file"
    #cat "${stdout_file}/1"

    ## TODO: Restore the variables (doesn't work currently because variables are printed using `env`)
    pash_redir_output echo "$$: (2) Recovering script variables from: $output_variable_file"
    source "$RUNTIME_DIR/pash_restore_fds.sh" "${output_variable_file}.fds" "${stdout_file}"
    source "$RUNTIME_DIR/pash_source_declare_vars.sh" "$output_variable_file"
    
elif [[ "$daemon_response" == *"UNSAFE:"* ]]; then
    pash_redir_output echo "$$: (2) Scheduler responded: $daemon_response"
    pash_redir_output echo "$$: (2) Executing command: $pash_speculative_command_id"
    ## Execute the command.
    ## KK 2023-06-01 Does `eval` work in general? We need to be precise
    ##               about which commands are unsafe to determine how to execute them.
    cmd="$(cat "$PASH_SPEC_NODE_DIRECTORY/$pash_speculative_command_id")"
    ## KK 2023-06-01 Not sure if this shellcheck warning must be resolved:
    ## > note: Double quote to prevent globbing and word splitting.
    # shellcheck disable=SC2086
    cmd_exit_code=$?
    if [ -z "${PASH_OLD_IFS+x}" ]; then
	unset IFS
    else
	IFS="$PASH_OLD_IFS"
    fi
    eval "$cmd"
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
unset cmd_exit_code
unset output_variable_file
unset cmd
unset stdout_file

## TODO: Also need to use wrap_vars maybe to `set` properly etc
