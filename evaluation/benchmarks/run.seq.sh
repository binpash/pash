#!/bin/bash

# FIXME: skip running if output file exists (using tee?)

## FIX: We should not have a set -e in a script that is supposed to be sourced.
# set -e

# time: print real in seconds, to simplify parsing
TIMEFORMAT="%3R" # %3U %3S"

if [[ -z "$PASH_TOP" ]]; then
  echo "Need to provide PASH_TOP, possibly $(git rev-parse --show-toplevel)" 1>&2
  exit 1
fi

oneliners(){
  seq_times_file="seq.res"
  seq_outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "oneliners/$seq_times_file" ]; then
    echo "skipping oneliners/$seq_times_file"
    return 0
  fi
  
  cd oneliners/

  cd ./input/
  ./setup.sh --full
  cd ..

  mkdir -p "$outputs_dir"

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

  touch "$seq_times_file"
  echo executing one-liners $(date) | tee -a "$seq_times_file"
  echo '' >> "$seq_times_file"

  for script_input in ${scripts_inputs[@]}
  do
    IFS=";" read -r -a script_input_parsed <<< "${script_input}"
    script="${script_input_parsed[0]}"
    input="${script_input_parsed[1]}"
    export IN="$PASH_TOP/evaluation/benchmarks/oneliners/input/$input"
    printf -v pad %30s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:30}

    seq_outputs_file="${outputs_dir}/${script}.${seq_outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$seq_outputs_file"; } 2>&1) | tee -a "$seq_times_file"
  done

  cd ..
}

unix50(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "unix50/${times_file}" ]; then
    echo "skipping unix50/${times_file}"
    return 0
  fi

  cd unix50/

  cd input/
  ./setup.sh
  cd ..

  mkdir -p "$outputs_dir"

  touch "$times_file"
  echo executing Unix50 $(date) | tee -a "$times_file"
  echo '' >> "$times_file"

  # FIXME this is the input prefix; do we want all to be IN 
  export IN_PRE=$PASH_TOP/evaluation/benchmarks/unix50/input

  for number in `seq 36`
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$times_file"
  done  
  cd ..
}

web-index(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "web-index/${times_file}" ]; then
    echo "skipping web-index/${times_file}"
    return 0
  fi

  cd web-index/

  cd input/
  ./setup.sh
  cd ..

  mkdir -p "$outputs_dir"
  
  touch "$times_file"
  echo executing web index $(date) | tee -a "$times_file"
  export IN=$PASH_TOP/evaluation/benchmarks/web-index/input/1000.txt
  export WEB_INDEX_DIR=$PASH_TOP/evaluation/benchmarks/web-index/input
  export WIKI=$PASH_TOP/evaluation/benchmarks/web-index/input/
  outputs_file="${outputs_dir}/web-index.${outputs_suffix}"
  echo web-index.sh: $({ time ./web-index.sh > "${outputs_file}"; } 2>&1) | tee -a "$times_file"
  cd ..
}

max-temp(){
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  cd max-temp/

  echo '' > seq.res
  echo executing max temp $(date) | tee -a ./seq.res
  echo mex-temp.sh: $({ time ./max-temp.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  cd ..
}

analytics-mts(){
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  cd analytics-mts/

  cd input/
  ./setup.sh
  cd ..

  echo '' > seq.res
  echo executing MTS analytics $(date) | tee -a ./seq.res
  echo 1.sh: $({ time ./1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 2.sh: $({ time ./2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 3.sh: $({ time ./3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 4.sh: $({ time ./4.sh > /dev/null; } 2>&1) | tee -a ./seq.res
#FIXME: echo 5.sh: $({ time ./5.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  cd ..
}

poets(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "poets/${times_file}" ]; then
    echo "skipping poets/${times_file}"
    return 0
  fi

  cd poets/

  cd input/
  ./setup.sh
  cd ..

  mkdir -p "$outputs_dir"

  names_scripts=(
    "1syllable_words;6_4"
    "2syllable_words;6_5"
    "4letter_words;6_2"
    "bigrams_appear_twice;8.2_2"
    "bigrams;4_3"
    "compare_exodus_genesis;8.3_3"
    "count_consonant_seq;7_2"
    # "count_morphs;7_1"
    "count_trigrams;4_3b"
    "count_vowel_seq;2_2"
    "count_words;1_1"
    "find_anagrams;8.3_2"
    "merge_upper;2_1"
    "sort;3_1"
    "sort_words_by_folding;3_2"
    "sort_words_by_num_of_syllables;8_1"
    "sort_words_by_rhyming;3_3"
    # "trigram_rec;6_1"
    "uppercase_by_token;6_1_1"
    "uppercase_by_type;6_1_2"
    "verses_2om_3om_2instances;6_7"
    "vowel_sequencies_gr_1K;8.2_1"
    "words_no_vowels;6_3"
  )

  touch "$times_file"
  echo executing Unix-for-poets $(date) | tee -a "$times_file"
  echo '' >> "$times_file"

  for name_script in ${names_scripts[@]}
  do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"
    export IN="$PASH_TOP/evaluation/benchmarks/poets/input/genesis"
    printf -v pad %30s
    padded_script="${name}.sh:${pad}"
    padded_script=${padded_script:0:30}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$times_file"
  done
  cd ..
}


aliases(){
  echo executing aliases
  cd aliases/select
  echo tomp3: $({ time ./1.tomp3.sh > /dev/null; } 2>&1)
  echo unrtf: $({ time ./2.unrtf.sh > /dev/null; } 2>&1)

}
