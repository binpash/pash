#!/bin/bash

OUTPUT=${1:="pash.tar.gz"}
BRANCH=${2:="master"}

# TODO: Make the temp_repo_dir be variable and random named

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
tar -czf ../$OUTPUT pash
cd ../
rm -rf temp_repo_dir
