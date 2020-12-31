#!/bin/bash

set -ex

echo $(pwd)
REV=0

REV=$(git rev-parse --short HEAD)
cd ../../

git clone --depth 1 git@github.com:andromeda/pash.git

# TODO: for a clear release, remove all versioning artifacts
# cp -r pash release
# cd release
# rm -rf .gitignore .gitsubmodules .git
# cd ..

mv pash pash-shallow
tar -cvzf pash-shallow.tar.gz pash-shallow/ > /dev/null
# uncomment the following line to keep all versions
# mv pash.tar.gz get/pash-${REV}.tar.gz
# ln -sf ./pash-${REV}.tar.gz get/latest 
mv pash-shallow.tar.gz get/
ln -sf ./pash-shallow.tar.gz get/latest 
rm -rf pash-shallow

# in the future, we might want to have versions
# ln -s pash.tar.gz latest
