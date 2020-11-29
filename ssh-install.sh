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

echo "Preparing repo archive..."
mkdir temp_repo_dir
cd temp_repo_dir
git clone git@github.com:andromeda/pash.git
cd pash
git checkout $BRANCH
## This needs to happen since libdash repo is linked using ssh
## TODO: Change the libdash submodule url to https to be cloneable without ssh keys
git submodule init
git submodule update
cd ../
tar -czf pash.tar.gz pash

echo "Sending repo archive..."
scp -o StrictHostKeyChecking=no -i $PRIVATE_KEY pash.tar.gz $USER@$HOSTNAME:/home/$USER

cd ../
rm -rf temp_repo_dir

ssh -o StrictHostKeyChecking=no -i $PRIVATE_KEY $USER@$HOSTNAME /bin/bash <<'EOF'
tar -xzf pash.tar.gz
cd pash
./install.sh -p
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
cd compiler
./test_evaluation_scripts.sh
EOF

