#!/bin/bash
  
# Specify the folder where the .out files are located
folder="$DISH_TOP/evaluation/distr_benchmarks/oneliners/outputs"

# Loop through the files in the folder
num_workers=3
for script_distr_out in "$folder"/*distr.out; do
    # Extract the script name without the extension
    script_name=$(basename "$script_distr_out" .distr.out)
    for ((i = 1; i <= num_workers; i++)); do
        # get the corresponding .faults.$crashed_worker.out file
        crashed_worker="worker$i"
        script_faults_out="$folder/$script_name.faults_$crashed_worker.out"

        # Perform a diff between the two files
        echo "Comparing faults_$crashed_worker.out and distr.out for script $script_name.sh"
        if diff -q "$script_faults_out" "$script_distr_out"; then
            echo "Outputs are identical"
        else
            echo "Files are different. Differences are as follows:"
            diff -y "$script_faults_out" "$script_distr_out"
        fi
        echo "-------------------------------------------"
    done
    
done