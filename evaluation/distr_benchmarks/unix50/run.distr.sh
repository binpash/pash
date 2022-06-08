PASH_FLAGS='--width 8 --r_split'
export TIMEFORMAT=%R

if [[ "$1" == "--extended" ]]; then
    echo "Using extended input"
    export IN_PRE=/unix50/extended_input
  else
    export IN_PRE=/unix50
fi

unix50_bash(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"

  mkdir -p "$outputs_dir"

  touch "$times_file"
  cat "$times_file" >> "$times_file".d
  echo executing Unix50 $(date) | tee "$times_file"
  echo '' >> "$times_file"

  for number in `seq 36`
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$times_file"
  done
}


unix50_pash(){
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
  cat "$times_file" >> "$times_file".d
  echo executing Unix50 $(date) | tee "$times_file"
  echo '' >> "$times_file"

  for number in `seq 36`
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"
    single_time_file="${outputs_dir}/${script}.${time_suffix}"

    echo -n "${padded_script}" | tee -a "$times_file"
    { time "$PASH_TOP/pa.sh" $flags --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done  
}

unix50_bash

unix50_pash "$PASH_FLAGS" "par"

unix50_pash "$PASH_FLAGS --distributed_exec" "distr"
