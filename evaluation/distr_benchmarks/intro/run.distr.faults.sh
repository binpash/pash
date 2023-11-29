PASH_FLAGS='--width 8 --r_split'
export TIMEFORMAT=%R
export dict="$PASH_TOP/evaluation/distr_benchmarks/oneliners/input/dict.txt"
curl -sf 'http://ndr.md/data/dummy/dict.txt' | sort > $dict


intro_pash(){
  flags=${1:-$PASH_FLAGS}
  prefix=${2:-par}
  prefix=$prefix

  times_file="$prefix.res"
  outputs_suffix="$prefix.out"
  time_suffix="$prefix.time"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs_$prefix"

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  touch "$times_file"
  cat $times_file >> $times_file.d
  echo executing one-liners with $prefix pash with data $(date) | tee "$times_file"
  echo '' >> "$times_file"


  script="demo-spell"


  printf -v pad %30s
  padded_script="${script}.sh:${pad}"
  padded_script=${padded_script:0:30}

  outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
  pash_log="${pash_logs_dir}/${script}.pash.log"
  single_time_file="${outputs_dir}/${script}.${time_suffix}"

  echo -n "${padded_script}" | tee -a "$times_file"
  { time "$PASH_TOP/pa.sh" $flags --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
  cat "${single_time_file}" | tee -a "$times_file"
  
}

intro_faults() {
  # For faults, mock crash for all workers
  num_workers=3
  # it's important to set the timeout long enough for now to avoid the "crashed" worker coming back alive while its replacement does work
  # until it's fully supported! 
  timeout=100

  for ((i = 1; i <= num_workers; i++)); do
    crashed_worker="worker$i"
    echo Mocking crash for $crashed_worker with timeout of $timeout seconds
    echo ----------------------------------------------------------------
    intro_pash "$PASH_FLAGS --distributed_exec --worker_timeout 100 --worker_timeout_choice worker$i" "faults_$crashed_worker"
    # echo "Iteration $i"
    # Your loop body here
  done
}

outputs_dir="outputs"
rm -rf "$outputs"

intro_pash "$PASH_FLAGS --distributed_exec" "distr"

intro_faults
