#!/usr/bin/env bash

## Steps:
## - clone the repository to a temp directory (requires access to clone the repo)
## - Compress the repository
## - scp copy the compressed repository to the server (requires ssh keys to server)
## - ssh to the server (requires ssh keys to server)
##   + untar the compressed repository
##   + run install script (requires sudo access without password)
##   + cd compiler; ./test_evaluation_scripts.sh

set -e

if [ $# -lt 3 ]; then
 echo "Not enough arguments!"
 exit 1
fi

## The directory to save the results in
LOCAL_DIR=$1

## The path of the private key for ssh authentication
PRIVATE_KEY=$2

## The user and IP address or hostname of the server
HOSTNAME=$3

USER="ubuntu"
RESULT_DIR="/home/${USER}/pash/evaluation/results/eurosys_small"

mkdir -p "$LOCAL_DIR"
scp -o StrictHostKeyChecking=no -o 'ConnectionAttempts 10' -i $PRIVATE_KEY -r "${USER}@${HOSTNAME}:${RESULT_DIR}" "${LOCAL_DIR}"
