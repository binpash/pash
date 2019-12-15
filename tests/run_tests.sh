#!/bin/bash

#Declare a string array
TestsArray=("minimal_sort_2" \
            "minimal_grep_2")

for test_name in ${TestsArray[*]}; do
     ./diff_test.sh $test_name /tmp/distr_output 1 100
done
