#!/bin/bash

## This sources variables that were produced from `declare -p`

## TODO: Fix this to not source read only variables
## TODO: Does this work with arrays

## TODO: Fix this to not source pash variables so as to not invalidate PaSh progress

## TODO: Fix this filtering

filter_vars_file()
{
    cat "$1" | grep -v "^declare -\([A-Za-z]\|-\)* \(pash\|BASH\|LINENO\|EUID\|GROUPS\)" 
    # The extension below is done for the speculative pash
    # | grep -v "LS_COLORS"
}

## TODO: Error handling if the argument is empty?
if [ "$PASH_DEBUG_LEVEL" -eq 0 ]; then
        > /dev/null 2>&1 "$@"
else
    if [ "$PASH_REDIR" == '&2' ]; then
        >&2 source <(filter_vars_file "$1")
    else
        >>"$PASH_REDIR" 2>&1 source <(filter_vars_file "$1")
    fi
fi
