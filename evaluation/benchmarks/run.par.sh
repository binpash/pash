#!/bin/bash

# time: print real in seconds, to simplify parsing
TIMEFORMAT="%3R" # %3U %3S"

if [[ -z "$PASH_TOP" ]]; then
  echo "Need to provide PASH_TOP, possibly $(git rev-parse --show-toplevel)" 1>&2
  exit 1
fi

oneliners_pash(){
  par_times_file="par.res"
  par_outputs_suffix="par.out"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs"
  width=16
  if [ -e "oneliners/$par_times_file" ]; then
    echo "skipping oneliners/$par_times_file"
    return 0
  fi
  
  cd oneliners/

  cd ./input/
  ./setup.sh --full
  cd ..

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  scripts_inputs=(
    "nfa-regex;100M.txt"
    "sort;3G.txt"
    "top-n;1G.txt"
    "wf;3G.txt"
    "spell;1G.txt"
    "diff;3G.txt"
    "bi-grams;1G.txt"
    "set-diff;3G.txt"
    "sort-sort;1G.txt"
    "shortest-scripts;all_cmdsx100.txt"
  )

  touch "$par_times_file"
  echo executing one-liners with pash $(date) | tee -a "$par_times_file"
  echo '' >> "$par_times_file"

  for script_input in ${scripts_inputs[@]}
  do
    IFS=";" read -r -a script_input_parsed <<< "${script_input}"
    script="${script_input_parsed[0]}"
    input="${script_input_parsed[1]}"
    export IN="$PASH_TOP/evaluation/benchmarks/oneliners/input/$input"
    printf -v pad %30s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:30}

    par_outputs_file="${outputs_dir}/${script}.${par_outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"

    echo "${padded_script}" $({ time "$PASH_TOP/pa.sh" --r_split --dgsh_tee -d 1 -w "${width}" --log_file "${pash_log}" ${script}.sh > "$par_outputs_file"; } 2>&1) | tee -a "$par_times_file"
  done

  cd ..
}
