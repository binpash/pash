#!/usr/bin/env bash

##
## Works for bash
##

## Configuration:
## VARIABLE_FILE_PREFIX: Set this variable to determine a specific prefix for the variable file
## DEFAULT_SET_STATE: Set this variable to determine the safe set state "huB"
##

## Save the previous exit code
export PREVIOUS_SHELL_EC="$?"


## TODO: Remove this from here
VARIABLE_FILE_PREFIX=$PASH_TMP_PREFIX


## TODO: There is somethin wrong with this now
## Use BASH_ARGV to save the arguments 
declare -ax PREVIOUS_SHELL_ARGS
# echo "${BASH_ARGV[*]}"
curr_argc=${BASH_ARGC[0]}
prev_argc=${BASH_ARGC[1]}
# echo "$curr_argc $prev_argc"
reversed_shell_args=("${BASH_ARGV[@]:$curr_argc:$prev_argc}")
# echo "${reversed_shell_args[*]}"
for (( i=1 ; i<=$prev_argc ; i++ )); 
do
    rev_i=$(( $prev_argc - $i ))
    PREVIOUS_SHELL_ARGS[i]=${reversed_shell_args[$rev_i]}
done
# echo "${PREVIOUS_SHELL_ARGS[*]}"
export PREVIOUS_SHELL_ARGS
# declare 

## TODO: Save the set state
export PREVIOUS_SET_STATUS=$-

## Save the shell variables to a file
export PREVIOUS_VARIABLES_FILE="${VARIABLE_FILE_PREFIX:-/tmp}/variables_$RANDOM$RANDOM$RANDOM"
## TODO: Rename declare_vars to something else and maybe remove RUNTIME_DIR
source "$RUNTIME_DIR/pash_declare_vars.sh" "$PREVIOUS_VARIABLES_FILE"

## TODO: Rename declare_vars to something else and maybe remove RUNTIME_DIR
# echo "$$: (1) Bash set state at start of execution: $pash_previous_set_status"
source "$RUNTIME_DIR/pash_set_from_to.sh" "$PREVIOUS_SET_STATUS" "${DEFAULT_SET_STATE:-huB}"
# echo "$$: (1) Set state reverted to PaSh-internal set state: $-"
