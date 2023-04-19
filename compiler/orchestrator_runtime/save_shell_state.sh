#!/usr/bin/env bash

##
## Works for bash
##

## Configuration:
## DEFAULT_SET_STATE: Set this variable to determine the safe set state "huB"
##

## Save the previous exit code
export PREVIOUS_SHELL_EC="$?"


# ## TODO: There is somethin wrong with this now, it messes up with test16
# ## Use BASH_ARGV to save the arguments 
# declare -ax PREVIOUS_SHELL_ARGS
# # echo "${BASH_ARGV[*]}"
# curr_argc=${BASH_ARGC[0]}
# prev_argc=${BASH_ARGC[1]}
# # echo "$curr_argc $prev_argc"
# reversed_shell_args=("${BASH_ARGV[@]:$curr_argc:$prev_argc}")
# # echo "${reversed_shell_args[*]}"
# for (( i=1 ; i<=$prev_argc ; i++ )); 
# do
#     rev_i=$(( $prev_argc - $i ))
#     PREVIOUS_SHELL_ARGS[i]=${reversed_shell_args[$rev_i]}
# done
# # echo "${PREVIOUS_SHELL_ARGS[*]}"
# export PREVIOUS_SHELL_ARGS
# # declare 
export PREVIOUS_SHELL_ARGS=( "$@" )


## Store the current `set` status 
export PREVIOUS_SET_STATUS=$-
source "$RUNTIME_DIR/pash_set_from_to.sh" "$PREVIOUS_SET_STATUS" "${DEFAULT_SET_STATE:-huB}"
