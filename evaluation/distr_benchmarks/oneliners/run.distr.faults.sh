PASH_FLAGS='--width 8 --r_split'
export TIMEFORMAT=%R
export dict="$PASH_TOP/evaluation/distr_benchmarks/oneliners/input/dict.txt"
curl -sf 'http://ndr.md/data/dummy/dict.txt' | sort > $dict


scripts_inputs=(
      "nfa-regex;1G.txt"
      "sort;3G.txt"
      "top-n;3G.txt"
      "wf;3G.txt"
      "spell;3G.txt"
      "diff;3G.txt"
      "bi-grams;3G.txt"
      "set-diff;3G.txt"
      "sort-sort;3G.txt"
      "shortest-scripts;all_cmdsx100.txt"
  )

# scripts_num_subgraphs=(
#     "nfa-regex;1"
#     "sort;1"
#     "top-n;1"
#     "wf;1"
#     "spell;1"
#     "diff;2"
#     "bi-grams;1"
#     "set-diff;2"
#     "sort-sort;1"
#     "shortest-scripts;1"
# )
# declare -A num_subgraphs_map

# # Populate the associative array
# for num_subgraph in "${scripts_num_subgraphs[@]}"; do
#   IFS=";" read -r -a subgraph_info <<< "$num_subgraph"
#   script_name="${subgraph_info[0]}"
#   num_subgraphs="${subgraph_info[1]}"
#   num_subgraphs_map["$script_name"]=$num_subgraphs
# done

oneliners_pash(){
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

  for script_input in ${scripts_inputs[@]}
  do
    IFS=";" read -r -a script_input_parsed <<< "${script_input}"
    script="${script_input_parsed[0]}"
    input="${script_input_parsed[1]}"

    export IN="/oneliners/$input"
    export dict=

    printf -v pad %30s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:30}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"
    single_time_file="${outputs_dir}/${script}.${time_suffix}"

    echo -n "${padded_script}" | tee -a "$times_file"
    { time "$PASH_TOP/pa.sh" $flags --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done
}

oneliners_faults() {
  # For faults, mock crash for all workers
  num_workers=3
  # it's important to set the timeout long enough for now to avoid the "crashed" worker coming back alive while its replacement does work
  # until it's fully supported! 
  timeout=100

  for ((i = 1; i <= num_workers; i++)); do
    crashed_worker="worker$i"
    echo Mocking crash for $crashed_worker with timeout of $timeout seconds
    echo ----------------------------------------------------------------
    oneliners_pash "$PASH_FLAGS --distributed_exec --worker_timeout 100 --worker_timeout_choice worker$i" "faults_$crashed_worker"
    # echo "Iteration $i"
    # Your loop body here
  done
}

outputs_dir="outputs"
rm -rf "$outputs"

oneliners_pash "$PASH_FLAGS --distributed_exec" "distr"

oneliners_faults
