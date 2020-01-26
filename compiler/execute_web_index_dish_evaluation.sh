#!/bin/bash

eval_dir="../evaluation/"
directory="${eval_dir}/scripts/web-index/"
p1="${directory}/p1.sh"
p2="${directory}/p2.sh"

output_dir="/tmp/web-index-output/"
p1_out="${output_dir}/p1.out"
p2_out="${output_dir}/p2.out"
p3_out="${output_dir}/p3.out"

intermediary_dir="${eval_dir}/intermediary/"
results_dir="${eval_dir}/results/"

input_dir="${HOME}/wikipedia/"

## Make the temporary output dir
mkdir -p $output_dir

cp "$input_dir/index_h_100.txt" $p1_out

{ time ( cat $p1_out | /bin/bash $p2 > $p2_out ) ; } 2> >(tee "${results_dir}/web-index_p2_seq.time" >&2)

echo "Sequential pipeline has been executed successfully."

## Split the intermediate output files for all scripts
split -n l/2 -d $p1_out ${p1_out}_2_
split -n l/10 -d $p1_out ${p1_out}_10_
split -n l/2 -d $p2_out ${p2_out}_2_
split -n l/10 -d $p2_out ${p2_out}_10_

echo "Intermediate files have been successfully produced."

## Setup the _env files of the experiments accordingly
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/web-index_p2_1_env.sh
echo "WIKI=${input_dir}" >> ${intermediary_dir}/web-index_p2_1_env.sh
echo "WEB_INDEX_DIR=${directory}" >> ${intermediary_dir}/web-index_p2_1_env.sh
cp ${intermediary_dir}/web-index_p2_1_env.sh ${intermediary_dir}/web-index_p2_2_env.sh
cp ${intermediary_dir}/web-index_p2_1_env.sh ${intermediary_dir}/web-index_p2_10_env.sh
cp ${intermediary_dir}/web-index_p2_1_funs.sh ${intermediary_dir}/web-index_p2_2_funs.sh
cp ${intermediary_dir}/web-index_p2_1_funs.sh ${intermediary_dir}/web-index_p2_10_funs.sh

## Setup the fifo pipes needed in the example.

## Note: This takes about 7 seconds for 100! urls
# ./execute_compile_evaluation_script.sh "web-index_p2_1"
./execute_compile_evaluation_script.sh "web-index_p2_2"
./execute_compile_evaluation_script.sh "web-index_p2_10"

## TODO: Fill in the rest of p2 and p3 here

