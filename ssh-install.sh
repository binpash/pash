#!/bin/bash

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

BRANCH="master"
USER="ubuntu"

ARCHIVE="pash.tar.gz"
echo "Preparing repo archive: $ARCHIVE branch: $BRANCH..."
./scripts/clone_compress_repo.sh $ARCHIVE $BRANCH &> /tmp/clone_compress.log

echo "Sending repo archive..."
scp -o StrictHostKeyChecking=no -i $PRIVATE_KEY pash.tar.gz "${USER}@${HOSTNAME}:/home/${USER}"

ssh -o StrictHostKeyChecking=no -i $PRIVATE_KEY "${USER}@${HOSTNAME}" /bin/bash <<'EOF'
tar -xzf pash.tar.gz
cd pash
echo "At branch: $(git rev-parse --abbrev-ref HEAD)"
source ./install.sh -p
cd compiler
./test_evaluation_scripts.sh
EOF

