#!/bin/bash

unix50_dir="../evaluation/unix50/"
unix50_intermediary="${unix50_dir}/intermediary/"
intermediary_dir="../evaluation/intermediary/"
results_subdir="unix50"

rm -r $unix50_intermediary
mkdir -p $unix50_intermediary
mkdir -p $intermediary_dir
mkdir -p "../evaluation/results/${results_subdir}/"

## Make inputs larger and generate scripts and their envs
input_size_increase=10
python3 generate_unix50_scripts.py $unix50_dir $unix50_intermediary $input_size_increase

n_inputs=(
    # 1
    2
    # 4
    # 8
    # 16
    # 32
    # 64
)

## Note: Maybe we need to tune the configuration (fan-out, batch-size)
##       before some specific benchmarks
for unix50_pipeline in $(ls ${unix50_intermediary} | grep -v "_env" | cut -f 1 -d '.' | sort); do
    echo $unix50_pipeline
    exec_seq="-s"
    for n_in in "${n_inputs[@]}"; do

        ## Generate the intermediary script
        python3 generate_microbenchmark_intermediary_scripts.py \
                $unix50_intermediary $unix50_pipeline $n_in $intermediary_dir

        ## Execute the intermediary script
        ./execute_compile_evaluation_script.sh $exec_seq -e "${unix50_pipeline}" "${n_in}" "${results_subdir}"  > /dev/null
        # Only execute the sequential once
        exec_seq=""

    done
done
