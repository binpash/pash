#!/usr/bin/env bash

##
## Works for bash
##

## Configuration:
## DEFAULT_SET_STATE: Set this variable to determine the safe set state "huB"
##

##
## Necessary for bash:
## - Last exit code $?
## - set state $-
##


## Save the previous exit code
export PREVIOUS_SHELL_EC="$?"

## Store the current `set` status 
export PREVIOUS_SET_STATUS=$-
source "$RUNTIME_DIR/pash_set_from_to.sh" "$PREVIOUS_SET_STATUS" "${DEFAULT_SET_STATE:-huB}"
