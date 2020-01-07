#!/bin/bash

directory="../evaluation/scripts/max-temp/"
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

intermediary_dir="../evaluation/intermediary/"

## Make the temporary output dir
mkdir -p $output_dir

# TODO: Maybe time p1, p2?
/bin/bash $p1 > $p1_out
cat $p1_out | /bin/bash $p2 > $p2_out

## TODO: p3 is currently only working for 2005
cat $p2_out | /bin/bash $p3 > $p3_out
cat $p3_out | /bin/bash $p4 > $p4_out
cat $p4_out | /bin/bash $p5 > $p5_out

echo "Sequential pipeline has been executed successfully."

## Split the intermediate output files for all scripts
split -n 2 -d $p3_out ${p3_out}_2_
split -n 10 -d $p3_out ${p3_out}_10_
split -n 2 -d $p4_out ${p4_out}_2_
split -n 10 -d $p4_out ${p4_out}_10_

echo "Intermediate files have been successfully produced."

## Setup the _env files of the experiments accordingly
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_1_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_2_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_10_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_1_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_2_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_10_env.sh

## TODO: Rename and remake the scripts in distr/seq to p4, p5
## TODO: Rename the input file names in the scripts to match those 

# ./execute_compile_evaluation_script.sh "p4_1"
# ./execute_compile_evaluation_script.sh "p4_2"
# ./execute_compile_evaluation_script.sh "p4_10"
# ./execute_compile_evaluation_script.sh "p5_1"
# ./execute_compile_evaluation_script.sh "p5_2"
# ./execute_compile_evaluation_script.sh "p5_10"


## (Maybe) TODO: Create the intermediary seq, distr scripts

