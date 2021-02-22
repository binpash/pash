#!/usr/bin/env bash

## Make this a first time only script

FORK="andromeda"
BRANCH="main"

rm -rf pash
git clone https://github.com/$FORK/pash.git
cd pash
git checkout $BRANCH
echo "At branch: \$(git rev-parse --abbrev-ref HEAD) of \$(git remote get-url origin)"
sed -i 's#git@github.com:angelhof/libdash.git#https://github.com/angelhof/libdash/#g' .gitmodules
source scripts/install.sh -p
cd compiler
./test_evaluation_scripts.sh