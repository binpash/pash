#!/bin/bash

eval_dir="../evaluation/"
directory="${eval_dir}/scripts/web-index/"
p1="${directory}/p1.sh"
p2="${directory}/p2.sh"

output_dir="/tmp/web-index-output/"
p1_out="${output_dir}/p1.out"
p2_out="${output_dir}/p2.out"

intermediary_dir="${eval_dir}/intermediary/"
results_dir="${eval_dir}/results/"

## Make the temporary output dir
mkdir -p $output_dir

{ time ( /bin/bash $p1 > $p1_out ) ; } 2> >(tee "${results_dir}/p1_seq.time" >&2)

echo "Sequential pipeline has been executed successfully."

## Split the intermediate output files for all scripts
split -n l/2 -d $p1_out ${p1_out}_2_
split -n l/10 -d $p1_out ${p1_out}_10_

echo "Intermediate files have been successfully produced."

## Setup the _env files of the experiments accordingly
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p2_1_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p2_2_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p2_10_env.sh

## Setup the fifo pipes needed in the example.

rm shift{1,2,3} {1,2,3}grams
mkfifo shift{1,2,3} {1,2,3}grams

./execute_compile_evaluation_script.sh "p2_1"



