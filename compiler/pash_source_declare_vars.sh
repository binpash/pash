#!/bin/bash

## This sources variables that were produced from `declare -p`

## TODO: Fix this to not source read only variables
## TODO: Does this work with arrays

## TODO: Fix this to not source pash variables so as to not invalidate PaSh progress

## TODO: Fix this filtering

## TODO: Error handling if the argument is empty?
if [ "$PASH_REDIR" == '&2' ]; then
    >&2 source <(cat $1 | grep -v "BASH" | grep -v "LINENO" | grep -v "EUID" | grep -v "GROUPS" | grep -v "pash")
else
    >>"$PASH_REDIR" 2>&1 source <(cat $1 | grep -v "BASH" | grep -v "LINENO" | grep -v "EUID" | grep -v "GROUPS" | grep -v "pash")
fi
