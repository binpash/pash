#!/bin/bash

echo $(pwd)
REV=0

git pull
REV=$(git rev-parse --short HEAD)
cd ../../

# TODO: for a clear release, remove all versioning artifacts
# cp -r pash release
# cd release
# rm -rf .gitignore .gitsubmodules .git
# cd ..

tar -cvzf pash.tar.gz pash/
# uncomment the following line to keep all versions
# mv pash.tar.gz get/pash-${REV}.tar.gz
# ln -sf ./pash-${REV}.tar.gz get/latest 
mv pash.tar.gz get/
ln -sf ./pash.tar.gz get/latest 

# in the future, we might want to have versions
# ln -s pash.tar.gz latest
