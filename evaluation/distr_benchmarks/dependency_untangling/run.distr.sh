PASH_FLAGS='--width 6 --r_split'
export TIMEFORMAT=%R
export dict="$PASH_TOP/evaluation/benchmarks/oneliners/input/dict.txt"

names_scripts=(
    "MediaConv1;img_convert"
    "MediaConv2;to_mp3"
    "Program_Inference;proginf"
    "LogAnalysis1;nginx"
    "LogAnalysis2;pcap"
    # "Genomics_Computation;genomics"
    "AurPkg;pacaur"
    "FileEnc1;compress_files"
    "FileEnc2;encrypt_files"
  )

oneliners_bash() {
    seq_times_file="seq.res"
    seq_outputs_suffix="seq.out"
    outputs_dir="outputs"

    mkdir -p "$outputs_dir"

    touch "$seq_times_file"
    cat $seq_times_file > $seq_times_file.d
    echo executing one-liners $(date) | tee -a "$seq_times_file"
    echo '' > "$seq_times_file"

    for name_script in ${names_scripts[@]}
    do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"
    export IN=
    export OUT=

    printf -v pad %30s
    padded_script="${script}${pad}"
    padded_script=${padded_script:0:30}

    seq_outputs_file="${outputs_dir}/${script}.${seq_outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$seq_outputs_file"; } 2>&1) | tee -a "$seq_times_file"
    done
}

oneliners_pash(){
  flags=${1:-$PASH_FLAGS}
  prefix=${2:-par}

  times_file="$prefix.res"
  outputs_suffix="$prefix.out"
  time_suffix="$prefix.time"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs_$prefix"

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  touch "$times_file"
  cat $times_file > $times_file.d
  echo executing one-liners with $prefix pash $(date) | tee -a "$times_file"
  echo '' > "$times_file"

  for name_script in ${names_scripts[@]}
    do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"

    export IN=
    export OUT=
    
    printf -v pad %30s
    padded_script="${script}${pad}"
    padded_script=${padded_script:0:30}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"
    single_time_file="${outputs_dir}/${script}.${time_suffix}"

    echo -n "${padded_script}" | tee -a "$times_file"
    { time "$PASH_TOP/pa.sh" $flags --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done
}

# oneliners_bash
oneliners_pash "$PASH_FLAGS" "par"
# oneliners_pash "$PASH_FLAGS --distributed_exec" "distr"
