PASH_FLAGS='--width 8 --r_split'
export TIMEFORMAT=%R

export IN="/max-temp/temperatures.txt"

max-temp_bash(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
 
  mkdir -p "$outputs_dir"
  touch "$times_file"
  cat "$times_file" >> "$times_file".d
  echo executing max temp $(date) | tee "$times_file"
  outputs_file="${outputs_dir}/temp-analytics.${outputs_suffix}"
  echo "temp-analytics.sh: " $({ time ./temp-analytics.sh > "${outputs_file}"; } 2>&1) | tee -a "$times_file"
}

max-temp_pash(){
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
  echo executing max-temp with $prefix pash $(date) | tee "$times_file"
  echo '' >> "$times_file"

  outputs_file="${outputs_dir}/temp-analytics.${outputs_suffix}"
  pash_log="${pash_logs_dir}/temp-analytics.pash.log"
  single_time_file="${outputs_dir}/temp-analytics.${time_suffix}"
  
  echo -n "temp-analytics.sh:  " | tee -a "$times_file"
  { time "$PASH_TOP/pa.sh" $flags --log_file "${pash_log}" temp-analytics.sh > "$outputs_file"; } 2> "${single_time_file}"
  cat "${single_time_file}" | tee -a "$times_file"
}

max-temp_hadoop_streaming(){
    jarpath="/opt/hadoop-3.4.0/share/hadoop/tools/lib/hadoop-streaming-3.4.0.jar" # Adjust as required
    infile="$IN"
    output_dirs="/outputs/hadoop-streaming/max-temp"
    times_file="hadoopstreaming.res"

    hadoop fs -rm -r $output_dirs/max
    hadoop fs -rm -r $output_dirs/min
    hadoop fs -rm -r $output_dirs/average

    cd "hadoop-streaming/"
    echo executing max-temp with hadoopstreaming $(date) | tee "$times_file"
    echo '' >> "$times_file"
    echo "temp-analytics.sh: " $({ time . run_all.sh &> /dev/null; } 2>&1) | tee -a "$times_file"
    cd ".."
    mv "hadoop-streaming/$times_file" .
}

max-temp_bash

max-temp_pash "$PASH_FLAGS" "par_no_du"

max-temp_pash "$PASH_FLAGS --parallel_pipelines --parallel_pipelines_limit 24" "par"

max-temp_pash "$PASH_FLAGS --distributed_exec" "distr_no_du"

max-temp_pash "$PASH_FLAGS --parallel_pipelines --distributed_exec --parallel_pipelines_limit 24" "distr"

max-temp_hadoop_streaming
