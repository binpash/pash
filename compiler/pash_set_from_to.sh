#!/bin/bash

from_set=${1?From set not given}
to_set=${2?To set not given}

## Finds the difference of set variables (removing the c, s one since it cannot be actually set and unset)
pash_redir_output echo "From set: $from_set"
pash_redir_output echo "To set: $to_set"
IFS=',' read -r pash_set_to_remove pash_set_to_add <<<"$("$RUNTIME_LIBRARY_DIR/set-diff" "$from_set" "$to_set")"
pash_redir_output echo "To add: $pash_set_to_add"
pash_redir_output echo "To remove: $pash_set_to_remove"
pash_redir_all_output_always_execute set "-$pash_set_to_add"
pash_redir_all_output_always_execute set "+$pash_set_to_remove"
