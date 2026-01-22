#!/bin/bash

score=0
total=0

if [ -d output ]; then
    echo "output directory already exists, aborting"
    exit 1
fi

mkdir output

echo "LEXER/PARSER AUTOGRADER RESULTS"
echo

# check success cases
for i in right/*.lc; do
    file=$(basename $i)
    output=$(mktemp output/$file.XXXX)

    echo -n "$file: "

    ./Main $i >$output 2>&1
    if [ $? -eq 0 ]
    then
        let score+=1
        echo "1/1"
    else
        echo "0/1"
    fi
    
    let total+=1
done

# check failure cases
for i in wrong/*.lc; do
    file=$(basename $i)
    output=$(mktemp output/$file.XXXX)

    echo -n "$file: "

    ./Main $i >$output 2>&1
    if [ $? -eq 1 ]
    then
        let score+=1
        echo "1/1"
    else
        echo "0/1"
    fi

    let total+=1
done

echo
echo "TOTAL: $score / $total"
echo
echo "PROBLEM 1: XXX / 5"
echo
let total=total+5
echo "FINAL GRADE: $score + XXX / $total"
