#!/bin/bash

## This script is necessary to ensure that sourcing happens with bash
# source run.seq.sh
# source run.par.sh

rm oneliners/par.res oneliners/seq.res

compare_outputs(){
  dir=$1
  outputs=$(ls $dir | grep "seq" | sed 's/.seq.out$//')
  for out in $outputs;
  do
    seq_output="${dir}/${out}.seq.out"
    pash_output="${dir}/${out}.par.out"
    diff -q "$seq_output" "$pash_output"
  done
}

source "$PASH_TOP/scripts/utils.sh"

oneliners(){
  seq_times_file="seq.res"
  seq_outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "oneliners/$seq_times_file" ]; then
    echo "skipping oneliners/$seq_times_file"
    return 0
  fi
  
  cd oneliners/
  # we need to download the whole dataset to generate the small input as well
  install_deps_source_setup $1

  mkdir -p "$outputs_dir"

  scripts_inputs=(
      "nfa-regex;100M.txt"
      # "sort;3G.txt"
      # "top-n;1G.txt"
      # "wf;3G.txt"
      # "spell;1G.txt"
      # "diff;3G.txt"
      # "bi-grams;1G.txt"
      # "set-diff;3G.txt"
      # "sort-sort;1G.txt"
      # "shortest-scripts;all_cmdsx100.txt"
  )

  touch "$seq_times_file"
  echo executing one-liners $(date) | tee -a "$seq_times_file"
  echo '' >> "$seq_times_file"

  for script_input in ${scripts_inputs[@]}
  do
    IFS=";" read -r -a script_input_parsed <<< "${script_input}"
    script="${script_input_parsed[0]}"
    input="${script_input_parsed[1]}"
    # source the required variables from setup.sh
    source_var ${1:-"notsmall"} $input
    printf -v pad %30s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:30}

    seq_outputs_file="${outputs_dir}/${script}.${seq_outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$seq_outputs_file"; } 2>&1) | tee -a "$seq_times_file"
  done

  cd ..
}

oneliners_pash(){
  #times_file="par.res"
  outputs_suffix="par.out"
  time_suffix="par.time"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs"
  if [ -e "oneliners/$times_file" ]; then
    echo "skipping oneliners/$times_file"
    return 0
  fi
  
  cd oneliners/

  install_deps_source_setup $1

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  scripts_inputs=(
      "nfa-regex;1M.txt" # 100M.txt"
      # "sort;3G.txt"
      # "top-n;1G.txt"
      # "wf;3G.txt"
      # "spell;1G.txt"
      # "diff;3G.txt"
      # "bi-grams;1G.txt"
      # "set-diff;3G.txt"
      # "sort-sort;1G.txt"
      # "shortest-scripts;all_cmdsx100.txt"
  )

  width_inputs=(
    2
    8
  )

  flags_inputs=(
    ""
    "--bash"
  )

  touch "$times_file"
  echo executing one-liners with pash $(date) | tee -a "$times_file"
  echo '' >> "$times_file"

  width

  for script_input in ${scripts_inputs[@]}
  do
    for width in ${width_inputs[@]}
    do
      for flags in ${flags_inputs[@]}
      do
        times_file="$par${flags}${width}.res"
        IFS=";" read -r -a script_input_parsed <<< "${script_input}"
        script="${script_input_parsed[0]}"
        input="${script_input_parsed[1]}"
        source_var $1 $input
        printf -v pad %30s
        padded_script="${script}.sh:${pad}"
        padded_script=${padded_script:0:30}

        outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
        pash_log="${pash_logs_dir}/${script}.pash.log"
        single_time_file="${outputs_dir}/${script}.${time_suffix}"

        echo -n "${padded_script}" | tee -a "$times_file"
        { time "$PASH_TOP/pa.sh" -w "${width}" "${flags}" $PASH_FLAGS --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
        cat "${single_time_file}" | tee -a "$times_file"
      done
    done
  done

  cd ..
}

oneliners --full

compare_outputs "oneliners/outputs"