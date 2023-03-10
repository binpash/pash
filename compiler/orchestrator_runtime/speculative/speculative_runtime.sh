#!/bin/bash


## TODO: Ask the scheduler to let us know when a command has been committed and what is its exit code.
## TODO: Define the client in pash_spec_init_setup (which should be sourced by pash_init_setup)

## TODO: Then we need to extend the scheduler to also support this protocol (unix sockets only) and 
##       Respond when the command is actually done.

export pash_speculative_command_id=$1

pash_redir_output echo "$$: (2) Before asking the scheduler for cmd: ${pash_speculative_command_id} exit code..."
## Send and receive from daemon
msg="Wait:${pash_speculative_command_id}"
daemon_response=$(pash_spec_communicate_scheduler "$msg") # Blocking step, daemon will not send response until it's safe to continue

## Receive an exit code
if [[ "$daemon_response" == *"OK:"* ]]; then
    response_args=($daemon_response)
    cmd_exit_code=${response_args[1]}
elif [ -z "$daemon_response" ]; then
    ## Trouble... Daemon crashed, rip
    pash_redir_output echo "$$: ERROR: (2) Scheduler crashed!"
    exit 1
else
    pash_redir_output echo "$$: ERROR: (2) Scheduler responded garbage ${daemon_response}!"
    exit 1
fi


pash_redir_output echo "$$: (2) Scheduler returned exit code: ${cmd_exit_code} for cmd with id: ${pash_speculative_command_id}."


## TODO: Figure out if this exits properly (prob not)
pash_runtime_final_status=${cmd_exit_code}

## TODO: Also need to use wrap_vars maybe to `set` properly etc
