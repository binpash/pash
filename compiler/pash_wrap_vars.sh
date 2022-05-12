#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

script_source="$1"

#ONLY WAY OUT IS TO TREAT EXEC in special way

## Recover the `set` state of the previous shell
# pash_redir_output echo "$$: (3) Previous BaSh set state: $pash_previous_set_status"
# pash_redir_output echo "$$: (3) PaSh-internal set state of current shell: $-"
export pash_current_set_state=$-
source "$RUNTIME_DIR/pash_set_from_to.sh" "$pash_current_set_state" "$pash_previous_set_status"
pash_redir_output echo "$$: (3) Reverted to BaSh set state: $-"

## Recover the input arguments of the previous script
## Note: We don't need to care about wrap_vars arguments because we have stored all of them already.
#
# This variable stores arguments as a space-separated stirng, so we need to
# unquote it and to split it into multiple strings by shell's field splitting.
# shellcheck disable=SC2086
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
    source "${script_source}"
    internal_exec_status=$?
    ## Make sure that any input argument changes are propagated outside
    export pash_input_args="$@"
    (exit "$internal_exec_status")
}
else 
{
    source "${script_source}"
    internal_exec_status=$?
    ## Make sure that any input argument changes are propagated outside
    export pash_input_args="$@"
    (exit "$internal_exec_status")
}
fi
