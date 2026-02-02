#!/bin/bash

from_set=${1?From set not given}
to_set=${2?To set not given}

## Finds the difference of set variables (removing the c, s one since it cannot be actually set and unset)
IFS=',' read -r __jit_set_to_remove __jit_set_to_add <<<"$("$RUNTIME_LIBRARY_DIR/set-diff" "$from_set" "$to_set")"
pash_redir_output echo "Switching set (from,to,add,remove): ($from_set,$to_set,$__jit_set_to_add,$__jit_set_to_remove)"
pash_redir_all_output_always_execute set "-$__jit_set_to_add"
pash_redir_all_output_always_execute set "+$__jit_set_to_remove"
