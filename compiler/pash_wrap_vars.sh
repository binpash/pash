#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

input_vars_file=${1??Input var file not given}

## Recover the variables that the previous shell had access to
## WARNING: This has to happen first, so that the variables underneath overwrite it.
pash_redir_output echo "$$: (3) Recovering variables from: $input_vars_file"
source "$RUNTIME_DIR/pash_source_declare_vars.sh" $input_vars_file

## Set the output file variables
## WARNING: This has to happen afterwards to avoid 
output_vars_file=${2?Output var file not given}
output_set_file=${3?Output set file not given}
script_source="${@:4}"

## Recover the `set` state of the previous shell
# pash_redir_output echo "$$: (3) Previous BaSh set state: $pash_previous_set_status"
# pash_redir_output echo "$$: (3) PaSh-internal set state of current shell: $-"
pash_current_set_state=$-
source "$RUNTIME_DIR/pash_set_from_to.sh" "$pash_current_set_state" "$pash_previous_set_status"
pash_redir_output echo "$$: (3) Reverted to BaSh set state: $-"

## Recover the input arguments of the previous script
previous_args="$@"
set -- $pash_input_args
pash_redir_output echo "$$: (3) Reverted to BaSh input arguments: $@"

## Execute the script
pash_redir_output echo "$$: (4) Restoring previous exit code: ${pash_previous_exit_status}"
pash_redir_output echo "$$: (4) Will execute script in ${script_source}:"
pash_redir_output cat "${script_source}"
## Note: We run the `exit` in a checked position so that we don't simply exit when we are in `set -e`.
if (exit "$pash_previous_exit_status")
then 
{
    ## Old way of executing the script but it doesn't properly work because the position is unchecked and therefore `set -e` doesn't behave as expected.
    # source "${script_source}" && internal_exec_status=$? || internal_exec_status=$?
    source "${script_source}"
    internal_exec_status=$?
}
else 
{
    source "${script_source}"
    internal_exec_status=$?
}
fi
pash_exec_status=${internal_exec_status}
pash_redir_output echo "$$: (5) BaSh script exited with ec: $pash_exec_status"

## Make sure that any input argument changes are propagated outside
export pash_input_args="$@"

## Recover the previous args of the script
set -- $previous_args

## Save the current set options to a file so that they can be recovered
pash_final_set_vars=$-
pash_redir_output echo "$$: (5) Writing current BaSh set state to: $output_set_file"
pash_redir_output echo "$$: (5) Current BaSh shell: $-"
echo "$pash_final_set_vars" > "$output_set_file"

## Revert to the old set state to avoid spurious fails
source "$RUNTIME_DIR/pash_set_from_to.sh" "$-" "$pash_current_set_state"
pash_redir_output echo "$$: (5) Reverted to PaSh set state to: $-"


## Save the current variables
source "$RUNTIME_DIR/pash_declare_vars.sh" $output_vars_file
pash_redir_output echo "$$: (5) Exiting from BaSh with BaSh status: $pash_exec_status"
(exit "$pash_exec_status")
