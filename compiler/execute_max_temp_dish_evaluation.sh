#!/bin/bash

directory="../evaluation/scripts/max-temp/"
p1="${directory}/p1.sh"
p2="${directory}/p2.sh"
p3="${directory}/p3.sh"
p4="${directory}/p4.sh"

output_dir="/tmp/max-temp-output/"
p1_out="${output_dir}/p1.out"
p2_out="${output_dir}/p2.out"
p3_out="${output_dir}/p3.out"
p4_out="${output_dir}/p4.out"

## Make the temporary output dir
mkdir -p $output_dir

# TODO: Maybe time p1, p2?
/bin/bash $p1 > $p1_out
cat $p1_out | /bin/bash $p2 > $p2_out
cat $p2_out | /bin/bash $p3 > $p3_out
cat $p3_out | /bin/bash $p4 > $p4_out

# ## We assume that each evaluation script has a sequential, a
# ## distributed, and an environment
# experiment=$1

# directory="../evaluation/"
# results="${directory}results/"
# prefix="${directory}${experiment}"
# env_file="${prefix}_env.sh"
# seq_script="${prefix}_seq.sh"
# distr_script="${prefix}_distr.sh"

# echo "Environment:"
# # cat $env_file
# . $env_file
# export $(cut -d= -f1 $env_file)

# echo "Sequential:"
# cat $seq_script | grep -v "rm -f" | grep -v "mkfifo"
# { time /bin/bash $seq_script > /tmp/seq_output ; } 2> >(tee "${results}${experiment}_seq.time" >&2)

# echo "Distributed:"
# { time ./translate_script.sh $seq_script $distr_script ; } 2> >(tee "${results}${experiment}_compile_distr.time" >&2)





