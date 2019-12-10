#!/bin/bash

## The first argument is the script to execute
script=$1

## The second argument is whether to compute the output difference
comput_diff=$2

## The third argument is a possible file with environment variable mappings
env=$3

json="/tmp/temp_script.json"
output_script="/tmp/temp_distr_script.sh"
seq_output="/tmp/seq_output"
distr_output="/tmp/distr_output"

## First parse the script_to_json
../parser/parse_to_json.native $script > $json

echo "Script:"
cat $script

# echo "Sequential:"
# time /bin/bash $script > $seq_output

## TODO: This currently assumes that there is only one IR in the script
echo "Distributed:"
python3 dish.py $json &&
time python3 distr_plan.py /tmp/dish_temp_ir_file1 $output_script $distr_output
# time /bin/bash $output_script

if [ "$comput_diff" = true ] ; then
    echo "Computing output difference... This will take some time..."
    diff -q $seq_output $distr_output
    if [[ $? == "0" ]]
    then
        echo "Output files are the same :)"
    else
        echo "ERROR: Output files are not the same"
    fi
fi
