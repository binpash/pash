PASH_FLAGS='--width 8 --r_split'
export TIMEFORMAT=%R

names_scripts=(
    "MediaConv1;img_convert"
    "MediaConv2;to_mp3"
    # "Program_Inference;proginf"
    "LogAnalysis1;nginx"
    "LogAnalysis2;pcap"
    # "Genomics_Computation;genomics"
    # "AurPkg;pacaur"
    "FileEnc1;compress_files"
    "FileEnc2;encrypt_files"
  )

dependency_untangling_bash() {
  outputs_dir="outputs"
  times_file="seq.res"
  outputs_suffix="seq.out"

  rm -rf input/output
  mkdir -p "$outputs_dir"

  touch "$times_file"
  cat "$times_file" > "$times_file".d
  echo executing dependency_untangling $(date) | tee "$times_file"
  echo '' >> "$times_file"
  
  export IN= 
  for name_script in ${names_scripts[@]}
  do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"
    printf -v pad %30s
    padded_script="${name}.sh:${pad}"
    padded_script=${padded_script:0:30}
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$times_file"
  done
}


dependency_untangling_pash() {
  flags=${1:-$PASH_FLAGS}
  prefix=${2:-par}

  times_file="$prefix.res"
  outputs_suffix="$prefix.out"
  time_suffix="$prefix.time"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs_$prefix"

  rm -rf input/output/

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  touch "$times_file"
  cat "$times_file" > "$times_file".d
  echo executing dependency_untangling with pash $(date) | tee "$times_file"
  echo '' >> "$times_file"
  
  export IN=
  for name_script in ${names_scripts[@]}
  do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"
    printf -v pad %30s
    padded_script="${name}.sh:${pad}"
    padded_script=${padded_script:0:30}
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"
    single_time_file="${outputs_dir}/${script}.${time_suffix}"
    
    echo -n "${padded_script}" | tee -a "$times_file"
    { time "$PASH_TOP/pa.sh" $flags --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done
}

dependency_untangling_bash

dependency_untangling_pash "$PASH_FLAGS" "par_no_du"

dependency_untangling_pash "$PASH_FLAGS --parallel_pipelines --parallel_pipelines_limit 24" "par"

dependency_untangling_pash "$PASH_FLAGS --distributed_exec" "distr_no_du"

dependency_untangling_pash "$PASH_FLAGS --parallel_pipelines --distributed_exec --parallel_pipelines_limit 24" "distr"

