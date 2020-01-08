#!/bin/bash

eval_dir="../evaluation/"
directory="${eval_dir}/scripts/max-temp/"
p1="${directory}/p1.sh"
p2="${directory}/p2.sh"
p3="${directory}/p3.sh"
p4="${directory}/p4.sh"
p5="${directory}/p5.sh"

output_dir="/tmp/max-temp-output/"
p1_out="${output_dir}/p1.out"
p2_out="${output_dir}/p2.out"
p3_out="${output_dir}/p3.out"
p4_out="${output_dir}/p4.out"
p5_out="${output_dir}/p5.out"

intermediary_dir="${eval_dir}/intermediary/"
results_dir="${eval_dir}/results/"

## Make the temporary output dir
mkdir -p $output_dir

time { /bin/bash $p1 > $p1_out ; } 2> >(tee "${results_dir}/p1_seq.time" >&2)
time { cat $p1_out | /bin/bash $p2 > $p2_out ; } 2> >(tee "${results_dir}/p2_seq.time" >&2)

## TODO: p3 is currently only working for 2005
time { cat $p2_out | /bin/bash $p3 > $p3_out ; } 2> >(tee "${results_dir}/p3_seq.time" >&2)
time { cat $p3_out | /bin/bash $p4 > $p4_out ; } 2> >(tee "${results_dir}/p4_seq.time" >&2)
time { cat $p4_out | /bin/bash $p5 > $p5_out ; } 2> >(tee "${results_dir}/p5_seq.time" >&2)

echo "Sequential pipeline has been executed successfully."

## Split the intermediate output files for all scripts
split -n l/2 -d $p3_out ${p3_out}_2_
split -n l/10 -d $p3_out ${p3_out}_10_
split -n l/2 -d $p4_out ${p4_out}_2_
split -n l/10 -d $p4_out ${p4_out}_10_

echo "Intermediate files have been successfully produced."

## Setup the _env files of the experiments accordingly
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_1_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_2_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_10_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_1_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_2_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_10_env.sh

./execute_compile_evaluation_script.sh "p4_1"
./execute_compile_evaluation_script.sh "p4_2"
./execute_compile_evaluation_script.sh "p4_10"
./execute_compile_evaluation_script.sh "p5_1"
./execute_compile_evaluation_script.sh "p5_2"
./execute_compile_evaluation_script.sh "p5_10"
