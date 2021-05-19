#!/usr/bin/env bash
set -e
cd $(dirname "$0")
TIMEFORMAT="%3R" # %3U %3S"
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
cd $PASH_TOP/evaluation/benchmarks/analytics-mts
cd input/
bash setup.sh $1
cd ..
export IN=input/in.csv
rm -f *.res
rm -rf outputs
rm -rf pash_logs
mkdir -p outputs
# time files
seq_times_file="seq.res"
par_times_file="par.res"
times_nc_file="par.nc.res"
# suffixes for output
seq_outputs_suffix="seq.out"
outputs_suffix="par.out"
outputs_nc_suffix="par.nc.out"
# outputs
outputs_dir="outputs"
pash_logs_dir="pash_logs"
mkdir -p $pash_logs_dir
mkdir -p $outputs_dir
for i in {1..4}
do
    script=$i
    echo "Executing $script..."
    seq_outputs_file="${outputs_dir}/${script}.${seq_outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.${outputs_suffix}"
    pash_nc_log="${pash_logs_dir}/${script}.${outputs_nc_suffix}"
    # SEQUENTIAL 
    echo "Executing the script with bash"
    echo "${script}" $({ time bash ${script}.sh > "$seq_outputs_file"; } 2>&1) | tee -a "$seq_times_file"
    # NO CAT SPLIT
    outputs_file_nc="${outputs_dir}/${script}.${outputs_nc_suffix}"
    echo "Executing the script with pash -w 16 without the cat-split optimization"
    echo "${script}" $({ time $PASH_TOP/pa.sh -d 1 -w 16 --log_file ${pash_nc_log} --no_cat_split_vanish ${script}.sh > ${outputs_file_nc} ; } 2>&1) | tee -a "${times_nc_file}"
    diff -s $seq_outputs_file $outputs_file_nc | head
    # PASH FULL
    echo "Executing the script with pash -w 16"
    par_outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    echo "${script}" $({ time $PASH_TOP/pa.sh  -d 1 -w 16 --log_file ${pash_log} ${script}.sh > "$par_outputs_file"; } 2>&1) | tee -a     "$par_times_file"
    diff -s $par_outputs_file $outputs_file_nc | head
done
paste seq.res par.nc.res par.res
