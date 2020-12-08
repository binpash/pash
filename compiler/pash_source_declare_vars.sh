#!/bin/bash

## This sources variables that were produced from `declare -p`

## TODO: Fix this to not source read only variables
## TODO: Does this work with arrays

## TODO: Error handling if the argument is empty?
if [ "$PASH_REDIR" == '&2' ]; then
    >&2 source $1
else
    >>"$PASH_REDIR" 2>&1 source $1
fi
