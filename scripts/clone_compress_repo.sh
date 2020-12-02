#!/bin/bash

OUTPUT=${1:="pash.tar.gz"}
BRANCH=${2:="master"}

# TODO: Make the temp_repo_dir be variable and random named

mkdir temp_repo_dir
cd temp_repo_dir
git clone --recursive git@github.com:andromeda/pash.git
cd pash
git checkout $BRANCH
cd ../
tar -czf ../$OUTPUT pash
cd ../
rm -rf temp_repo_dir
