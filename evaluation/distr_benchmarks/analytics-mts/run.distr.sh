PASH_FLAGS='--width 8 --r_split'
export TIMEFORMAT=%R

if [[ "$1" == "--small" ]]; then
    export IN="/analytics-mts/in_small.csv"
else
    export IN="/analytics-mts/in.csv"
fi

analytics-mts_bash(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"

  mkdir -p "$outputs_dir"

  touch "$times_file"
  cat "$times_file" >> "$times_file".d
  echo executing MTS analytics $(date) | tee "$times_file"
  echo '' >> "$times_file"
  ## FIXME 5.sh is not working yet
  for number in `seq 4`
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}
    # select the respective input
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$times_file"
  done
}

analytics-mts_pash(){
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
  echo executing MTS analytics with pash $(date) | tee "$times_file"
  echo '' >> "$times_file"
  ## FIXME 5.sh is not working yet
  for number in `seq 4`
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"
    single_time_file="${outputs_dir}/${script}.${time_suffix}"

    echo -n "${padded_script}" | tee -a "$times_file"
    { time "$PASH_TOP/pa.sh" $flags   --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done
}

analytics-mts_hadoopstreaming(){
  jarpath="/opt/hadoop-3.4.0/share/hadoop/tools/lib/hadoop-streaming-3.4.0.jar" # Adjust as required
  times_file="hadoopstreaming.res"
  outputs_suffix="hadoopstreaming.out"
  outputs_dir="/outputs/hadoop-streaming/analytics-mts"

  cd "hadoop-streaming/"

  hdfs dfs -rm -r "$outputs_dir"
  hdfs dfs -mkdir -p "$outputs_dir"

  touch "$times_file"
  cat "$times_file" >> "$times_file".d
  echo executing analytics-mts $(date) | tee "$times_file"
  echo '' >> "$times_file"

  COUNTER=1
  while IFS= read -r line; do
      printf -v pad %20s
      padded_script="${COUNTER}.sh:${pad}"
      padded_script=${padded_script:0:20}

      echo "${padded_script}" $({ time { eval $line &> /dev/null; } } 2>&1) | tee -a "$times_file"
      COUNTER=$(( COUNTER + 1 ))
  done <"run_all.sh"
  cd ".."
  mv "hadoop-streaming/$times_file" .
}

analytics-mts_bash

analytics-mts_pash "$PASH_FLAGS" "par"

analytics-mts_pash "$PASH_FLAGS --distributed_exec" "distr"

analytics-mts_hadoopstreaming

