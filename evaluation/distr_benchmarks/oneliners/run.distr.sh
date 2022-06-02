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

oneliners_bash() {
    outputs_dir="outputs"
    rep=${1:-rep3}
    seq_times_file=$rep"_seq.res"
    seq_outputs_suffix=$rep"_seq.out"

    mkdir -p "$outputs_dir"

    touch "$seq_times_file"
    cat $seq_times_file > $seq_times_file.d
    echo executing one-liners $(date) | tee "$seq_times_file"
    echo '' >> "$seq_times_file"

    for script_input in ${scripts_inputs[@]}
    do
    IFS=";" read -r -a script_input_parsed <<< "${script_input}"
    script="${script_input_parsed[0]}"
    input="${script_input_parsed[1]}"

    export IN=/$rep\_$input
    export dict=

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
  rep=${3:-rep3}
  prefix=$prefix\_$rep

  times_file="$prefix.res"
  outputs_suffix="$prefix.out"
  time_suffix="$prefix.time"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs_$prefix"

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  touch "$times_file"
  cat $times_file > $times_file.d
  echo executing one-liners with $prefix pash with data $rep $(date) | tee "$times_file"
  echo '' >> "$times_file"

  for script_input in ${scripts_inputs[@]}
  do
    IFS=";" read -r -a script_input_parsed <<< "${script_input}"
    script="${script_input_parsed[0]}"
    input="${script_input_parsed[1]}"

    export IN=/$rep\_$input
    export dict=

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

# oneliners_bash "rep1"
oneliners_bash "rep3"

# oneliners_pash "$PASH_FLAGS" "par" "rep1"
oneliners_pash "$PASH_FLAGS" "par" "rep3"

# oneliners_pash "$PASH_FLAGS --distributed_exec" "distr" "rep1"
oneliners_pash "$PASH_FLAGS --distributed_exec" "distr" "rep3"
