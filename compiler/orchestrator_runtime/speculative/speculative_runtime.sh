#!/bin/bash


## TODO: Ask the scheduler to let us know when a command has been committed and what is its exit code.
## TODO: Define the client in pash_spec_init_setup (which should be sourced by pash_init_setup)

## TODO: Then we need to extend the scheduler to also support this protocol (unix sockets only) and 
##       Respond when the command is actually done.

export pash_speculative_command_id=$1

echo "STUB: This would call the scheduler for command with id: ${pash_speculative_command_id}"

## TODO: Set this based on what the scheduler returns
pash_runtime_final_status=$?

## TODO: Also need to use wrap_vars maybe to `set` properly etc
