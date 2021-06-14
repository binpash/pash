#!/usr/bin/env bash
set -e
cd $(dirname "$0")
TIMEFORMAT="%3R" # %3U %3S"
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
cd $PASH_TOP/evaluation/benchmarks/oneliners
# setup 
cd input
# remove old input
./setup.sh -c 
./setup.sh --full
cd ..

scripts_inputs=(
"nfa-regex;100M.txt"
"sort;3G.txt"
"top-n;1G.txt"
"wf;3G.txt"
"spell;1G.txt"
"diff;3G.txt"
"bi-grams;1G.txt"
"set-diff;3G.txt"
#"sort-sort;1G.txt"
"shortest-scripts;all_cmdsx100.txt"
)

# Cleanup
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
# Generate the output folders
mkdir -p $pash_logs_dir
mkdir -p $outputs_dir

for script_input in ${scripts_inputs[@]}
do
  IFS=";" read -r -a script_input_parsed <<< "${script_input}"
  script="${script_input_parsed[0]}"
  input="${script_input_parsed[1]}"
  export IN="input/$input"
  printf -v pad %30s
  padded_script="${script}.sh:${pad}"
  padded_script=${padded_script:0:30}

  seq_outputs_file="${outputs_dir}/${script}.${seq_outputs_suffix}"
  pash_log="${pash_logs_dir}/${script}.${outputs_suffix}"
  pash_nc_log="${pash_logs_dir}/${script}.${outputs_nc_suffix}"
  # SEQUENTIAL 
  echo "Executing the script with bash"
  echo "${padded_script}" $({ time bash ${script}.sh > "$seq_outputs_file"; } 2>&1) | tee -a "$seq_times_file"
  # NO CAT SPLIT
  outputs_file_nc="${outputs_dir}/${script}.${outputs_nc_suffix}"
  echo "Executing the script with pash -w 16 without the cat-split optimization"
  echo "${padded_script}" $({ time $PASH_TOP/pa.sh -d 1 -w 16 --log_file ${pash_nc_log} --no_cat_split_vanish ${script}.sh > ${outputs_file_nc} ; } 2>&1) | tee -a "${times_nc_file}"
  diff -s $seq_outputs_file $outputs_file_nc | head
  # PASH FULL
  echo "Executing the script with pash -w 16"
  par_outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
  echo "${padded_script}" $({ time $PASH_TOP/pa.sh  -d 1 -w 16 --log_file ${pash_log} ${script}.sh > "$par_outputs_file"; } 2>&1) | tee -a     "$par_times_file"
  diff -s $par_outputs_file $outputs_file_nc | head
done
paste seq.res par.nc.res par.res
