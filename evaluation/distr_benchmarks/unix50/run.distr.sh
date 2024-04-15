PASH_FLAGS='--width 8 --r_split'
export TIMEFORMAT=%R
names_scripts=(
    1 2 3 4 5 6 7 8 9 10
    11 12 13 14 15 16 17 18 19 20 
    21 23 24 25 26 28 29
    30 31 32 33 34 35 36
  )
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

  for number in ${names_scripts[@]}
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh 2> /dev/null > "$outputs_file"; } 2>&1) | tee -a "$times_file"
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

  for number in ${names_scripts[@]}
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"
    single_time_file="${outputs_dir}/${script}.${time_suffix}"

    echo -n "${padded_script}" | tee -a "$times_file"
    { time "$PASH_TOP/pa.sh" $flags --log_file "${pash_log}" ${script}.sh 2> /dev/null > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done  
}

unix50_hadoopstreaming(){
  jarpath="/opt/hadoop-3.4.0/share/hadoop/tools/lib/hadoop-streaming-3.4.0.jar" # Adjust as required
  basepath="/unix50in" # Adjust as required
  times_file="hadoopstreaming.res"
  outputs_suffix="hadoopstreaming.out"
  outputs_dir="/outputs/hadoop-streaming/unix50"

  cd "hadoop-streaming/"

  hdfs dfs -rm -r "$outputs_dir"
  hdfs dfs -mkdir -p "$outputs_dir"

  touch "$times_file"
  cat "$times_file" >> "$times_file".d
  echo executing Unix50 $(date) | tee "$times_file"
  echo '' >> "$times_file"

  COUNTER=0
  while IFS= read -r line; do
      COUNTER=$(( COUNTER + 1 ))
      if [[ ! " ${names_scripts[*]} " =~ " ${COUNTER} " ]]; then
        continue # skip if array doesn't have this value
      fi
      printf -v pad %20s
      padded_script="${COUNTER}.sh:${pad}"
      padded_script=${padded_script:0:20} 
      echo "${padded_script}" $({ time { eval $line &> /dev/null; } } 2>&1) | tee -a "$times_file"
  done <"run_all.sh"
  cd ".."
  mv "hadoop-streaming/$times_file" .
}

unix50_bash

unix50_pash "$PASH_FLAGS" "par"

unix50_pash "$PASH_FLAGS --distributed_exec" "distr"

unix50_hadoopstreaming
