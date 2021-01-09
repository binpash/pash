#!/bin/bash

## TODO: Measure separate speedup for preprocessing and computation

## TODO: Maybe also measure the whole pipeline together (if it works)

eval_dir="../evaluation/"
directory="${eval_dir}/scripts/max-temp/"
p1="${directory}/p1.sh"
p2="${directory}/p2.sh"
p3="${directory}/p3.sh"
p4="${directory}/p4.sh"
p5="${directory}/p5.sh"

output_dir="$HOME/max-temp-output/"
p1_out="${output_dir}/p1.out"
p2_out="${output_dir}/p2.out"
p3_out="${output_dir}/p3.out"
p4_out="${output_dir}/p4.out"
p5_out="${output_dir}/p5.out"

intermediary_dir="${eval_dir}/intermediary/"
results_dir="${eval_dir}/results/"
 
## Make the temporary output dir
mkdir -p $output_dir

{ time ( /bin/bash $p1 > $p1_out ) ; } 2> >(tee "${results_dir}/p1_seq.time" >&2)
{ time ( cat $p1_out | /bin/bash $p2 > $p2_out ) ; } 2> >(tee "${results_dir}/p2_seq.time" >&2)
{ time ( cat $p2_out | /bin/bash $p3 > $p3_out ) ; } 2> >(tee "${results_dir}/p3_seq.time" >&2)
{ time ( cat $p3_out | /bin/bash $p4 > $p4_out ) ; } 2> >(tee "${results_dir}/p4_seq.time" >&2)
# { time ( cat $p4_out | /bin/bash $p5 > $p5_out ) ; } 2> >(tee "${results_dir}/p5_seq.time" >&2)

echo "Sequential pipeline has been executed successfully."

## Split the intermediate output files for all scripts
# split -n l/2 -d $p3_out ${p3_out}_2_
# split -n l/10 -d $p3_out ${p3_out}_10_
split -n l/16 -d $p3_out ${p3_out}_16_
# split -n l/2 -d $p4_out ${p4_out}_2_
# split -n l/10 -d $p4_out ${p4_out}_10_
split -n l/16 -d $p4_out ${p4_out}_16_

echo "Intermediate files have been successfully produced."
 
## Setup the _env files of the experiments accordingly
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_1_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_2_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_10_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p4_16_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_1_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_2_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_10_env.sh
echo "IN_DIR=${output_dir}" > ${intermediary_dir}/p5_16_env.sh

# ./execute_compile_evaluation_script.sh -s -e "p4" "2"
# rm -f /tmp/eager*
# ./execute_compile_evaluation_script.sh -e "p4" "10"
# rm -f /tmp/eager*
./execute_compile_evaluation_script.sh -e "p4" "16"
rm -f /tmp/eager*
# ./execute_compile_evaluation_script.sh -s -e "p5" "2"
# rm -f /tmp/eager*
# ./execute_compile_evaluation_script.sh -e "p5" "10"
# rm -f /tmp/eager*
./execute_compile_evaluation_script.sh -e "p5" "16"
rm -f /tmp/eager*

## Executing the first 3 stages with pash
# ./execute_compile_evaluation_script.sh -s -a "max_temp_p123" "2"
# rm -f /tmp/eager*
# ./execute_compile_evaluation_script.sh -a "max_temp_p123" "10"
# rm -f /tmp/eager*
./execute_compile_evaluation_script.sh -a "max_temp_p123" "16"
rm -f /tmp/eager*
