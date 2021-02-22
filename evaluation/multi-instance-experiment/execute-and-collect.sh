#!/usr/bin/env bash

INSTANCE_ID="${1:?No instance-id given.}"
KEY_PATH="${2:?No key path given.}"

## The directory to save the results in
LOCAL_DIR="${3:?No local results dir}"

## The user and IP address or hostname of the server
HOSTNAME="${4:?No hostname}"

EVAL_DIR="$PASH_TOP/evaluation/multi-instance-experiment/"

## TODO: This only works if the user's name is ubuntu (and also is tied to the value in experiment-script.sh)
##       Is there any better way to refactor this so that this is not an issue?
COMPLETION_FILE="/home/ubuntu/pash/done.txt"

## Execute script
"${PASH_TOP}/evaluation/multi-instance-experiment/ssh-script.sh" "${KEY_PATH}" "${HOSTNAME}" < "${EVAL_DIR}/experiment-script.sh" 

## Wait for it to complete
"${PASH_TOP}/scripts/ssh-bg-poller.sh" "${COMPLETION_FILE}" "${KEY_PATH}" "${HOSTNAME}"

## Collect results
"${EVAL_DIR}/collect-results-script.sh" "${LOCAL_DIR}" "${KEY_PATH}" "${HOSTNAME}"
