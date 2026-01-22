#!/bin/bash

if [ "$#" != "1" ]; then
    echo "Usage: $0 [hwXX]"
    exit 1
fi

if [ ! -d "$1/grading" ]; then
    echo "Couldn't find grading directory (looked in $1/grading)"
    exit 2
fi

cd $1/grading

errors=""
for s in `ls`; do
    echo "GRADING $s"
    (cd $s; make)
    if [ "$?" != "0" ]; then
        errors+=" $s"
    fi
done

echo
echo "There were errors for the following students:${errors}"
echo ${errors} >"$1/grading/errors.log"
