#!/bin/bash

from_set=${1?From set not given}
to_set=${2?To set not given}

## Finds the difference of set variables (removing the c, s one since it cannot be actually set and unset)
##
## TODO: This is very inefficient! It can probably happen in a single C command.
pash_set_to_add="$(comm -23 <(echo $to_set | sed -E 's/(.)/\1\n/g' | sort) <(echo $from_set | sed -E 's/(.)/\1\n/g' | sort) | grep -v 'c' | grep -v 's' || true )"
pash_set_to_remove="$(comm -13 <(echo $to_set | sed -E 's/(.)/\1\n/g' | sort) <(echo $from_set | sed -E 's/(.)/\1\n/g' | sort) | grep -v 'c' | grep -v 's' || true )"
pash_redir_output echo "To add: $pash_set_to_add"
pash_redir_output echo "To remove: $pash_set_to_remove"
pash_redir_all_output_always_execute set "-$pash_set_to_add"
pash_redir_all_output_always_execute set "+$pash_set_to_remove"