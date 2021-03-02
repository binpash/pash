#!/usr/bin/env bash

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

eval_dir="$PASH_TOP/evaluation/buses/"
results_dir="${eval_dir}/results/"

mkdir -p $results_dir

for i in 1 2 3 4
do
    script="${eval_dir}/${i}.sh"
    echo "Executing $script..."

    seq_output=/tmp/seq_output
    pash_width_16_no_cat_split_output=/tmp/pash_16_no_cat_split_output
    pash_width_16_output=/tmp/pash_16_output

    seq_time="${results_dir}/${i}_2_seq.time"
    pash_width_16_no_cat_split_time="${results_dir}/${i}_16_distr_auto_split_fan_in_fan_out.time"
    pash_width_16_time="${results_dir}/${i}_16_distr_auto_split.time"

    echo "Executing the script with bash..."
    { time /bin/bash $script > $seq_output ; } 2> >(tee "${seq_time}" >&2)

    echo "Executing the script with pash -w 16 without the cat-split optimization (log in: /tmp/pash_16_log)"
    { time $PASH_TOP/pa.sh -w 16 -d 1 --log_file /tmp/pash_16_no_cat_split_log --no_cat_split_vanish --output_time $script ; } 1> "$pash_width_16_no_cat_split_output" 2> >(tee "${pash_width_16_no_cat_split_time}" >&2)
    echo "Checking for output equivalence..."
    diff -s $seq_output $pash_width_16_no_cat_split_output | head

    echo "Executing the script with pash -w 16 (log in: /tmp/pash_16_log)"
    { time $PASH_TOP/pa.sh -w 16 -d 1 --log_file /tmp/pash_16_log --output_time $script ; } 1> "$pash_width_16_output" 2> >(tee "${pash_width_16_time}" >&2)
    echo "Checking for output equivalence..."
    diff -s $seq_output $pash_width_16_output | head

done
