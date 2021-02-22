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

if [ $# -lt 2 ]; then
 echo "Not enough arguments!"
 exit 1
fi

## The path of the private key for ssh authentication
PRIVATE_KEY=$1

## The user and IP address or hostname of the server
HOSTNAME=$2

USER="ubuntu"

ssh -o StrictHostKeyChecking=no -o 'ConnectionAttempts 10' -i $PRIVATE_KEY "${USER}@${HOSTNAME}" /bin/bash < "$PASH_TOP/evaluation/multi-instance-experiment/experiment-script.sh"
