#!/bin/bash

echo $(pwd)

cd pash
git pull
cd ..

# TODO: for a clear release, remove all versioning artifacts
# cp -r pash release
# cd release
# rm -rf .gitignore .gitsubmodules .git
# cd ..

tar -cvzf pash.tar.gz pash/
# in the future, we might want to have versions
# ln -s pash.tar.gz latest
