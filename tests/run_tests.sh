#!/bin/bash

#Declare a string array
TestsArray=("minimal_sort_2" \
            "minimal_grep_2")

## TODO: Exit if any diff_test fails
for test_name in ${TestsArray[*]}; do
     ./diff_test.sh $test_name /tmp/distr_output 1 100
done
