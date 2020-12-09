#!/bin/bash

from_set=${1?From set not given}
to_set=${2?To set not given}

## Finds the difference of set variables
pash_set_to_add="$(comm -23 <(echo $to_set | sed -E 's/(.)/\1\n/g' | sort) <(echo $from_set | sed -E 's/(.)/\1\n/g' | sort))"
pash_set_to_remove="$(comm -13 <(echo $to_set | sed -E 's/(.)/\1\n/g' | sort) <(echo $from_set | sed -E 's/(.)/\1\n/g' | sort))"
# pash_redir_output echo "To add: $pash_set_to_add"
# pash_redir_output echo "To remove: $pash_set_to_remove"
pash_redir_all_output set "-$pash_set_to_add"
pash_redir_all_output set "+$pash_set_to_remove"