#!/bin/bash

usecase_dir="../evaluation/usecases/shellcheck/"
intermediary_dir="../evaluation/intermediary/"

echo "Deleting eager intermediate files..."
rm -f /tmp/eager*
mkdir -p $intermediary_dir

n_inputs=(
    2
    # 4
    # 8
    # 16
    # 32
    # 64
)

name=distrotest

# Execute the sequential script on the first run only
exec_seq=""
for n_in in "${n_inputs[@]}"; do

    ## Generate the intermediary script
    cp "${usecase_dir}/${name}.sh" "${intermediary_dir}/${name}_${n_in}_seq.sh"
    cp "${usecase_dir}/${name}_env.sh" "${intermediary_dir}/${name}_${n_in}_env.sh"
    cp "${usecase_dir}/${name}_funs.sh" "${intermediary_dir}/${name}_${n_in}_funs.sh"

    ## Execute the intermediary script with eager
    ./execute_compile_evaluation_script.sh $exec_seq -e "${name}" "${n_in}"
    rm -f /tmp/eager*

    # Only execute the sequential once
    exec_seq=""

    ## Execute the intermediary script without eager
    # ./execute_compile_evaluation_script.sh $exec_seq "${microbenchmark}" "${n_in}"

    # ## Execute the intermediary script with the naive eager
    # ./execute_compile_evaluation_script.sh $exec_seq -n "${microbenchmark}" "${n_in}"
    # rm -f /tmp/eager*
done
