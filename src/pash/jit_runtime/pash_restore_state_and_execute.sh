#!/bin/bash

## File directory
RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")

## INPUT: Expects SCRIPT_TO_EXECUTE to be set

#ONLY WAY OUT IS TO TREAT EXEC in special way

## Recover the `set` state of the previous shell
# pash_redir_output echo "$$: (3) Previous BaSh set state: $pash_previous_set_status"
# pash_redir_output echo "$$: (3) PaSh-internal set state of current shell: $-"
export pash_current_set_state=$-
source "$RUNTIME_DIR/pash_set_from_to.sh" "$pash_current_set_state" "$pash_previous_set_status"
pash_redir_output echo "$$: (3) Reverted to BaSh set state: $-"

## Execute the script
pash_redir_output echo "$$: (4) Restoring previous exit code: ${pash_previous_exit_status}"
pash_redir_output echo "$$: (4) Will execute script in ${SCRIPT_TO_EXECUTE}:"
pash_redir_output cat "${SCRIPT_TO_EXECUTE}"

## Note: We run the `exit` in a checked position so that we don't simply exit when we are in `set -e`.
if (exit "$pash_previous_exit_status")
then 
{
    ## This works w.r.t. arguments because source does not change them if there are no arguments
    ## being given.
    source "${SCRIPT_TO_EXECUTE}"
}
else 
{
    source "${SCRIPT_TO_EXECUTE}"
}
fi
