#!/bin/bash

##
## This currently performs (5), i.e., reverting bash state to get back to pash mode.
##

## TODO: Use that for (1) too

output_vars_file=${1?Output var file not given}
output_set_file=${2?Output set file not given}

pash_exec_status=${internal_exec_status}
pash_redir_output echo "$$: (5) BaSh script exited with ec: $pash_exec_status"

## Save the current set options to a file so that they can be recovered
pash_final_set_vars=$-
pash_redir_output echo "$$: (5) Writing current BaSh set state to: $output_set_file"
pash_redir_output echo "$$: (5) Current BaSh shell: $-"
echo "$pash_final_set_vars" > "$output_set_file"

## Revert to the old set state to avoid spurious fails
source "$RUNTIME_DIR/pash_set_from_to.sh" "$-" "$pash_current_set_state"
pash_redir_output echo "$$: (5) Reverted to PaSh set state to: $-"


## Save the current variables
source "$RUNTIME_DIR/pash_declare_vars.sh" "$output_vars_file"
# pash_redir_output echo "$$: (5) Exiting from BaSh with BaSh status: $pash_exec_status"
# (exit "$pash_exec_status")
