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

filter_vars_file "$1"