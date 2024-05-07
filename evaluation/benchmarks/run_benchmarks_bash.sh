#!/bin/bash
## This script is necessary to ensure that sourcing happens with bash
# source run.seq.sh
# source run.par.sh

rm -f oneliners/par.res oneliners/seq.res

compare_outputs(){
  dir=$1
  outputs=$(ls $dir | grep "seq" | sed 's/.seq.out$//')
  for out in $outputs;
  do
    seq_output="${dir}/${out}.seq.out"
    pash_output="${dir}/${out}.par.out"
    diff -q "$seq_output" "$pash_output"
  done
}

source "$PASH_TOP/scripts/utils.sh"

oneliners(){
  seq_times_file="seq.res"
  seq_outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "oneliners/$seq_times_file" ]; then
    echo "skipping oneliners/$seq_times_file"
    return 0
  fi
  
  cd oneliners/
  # we need to download the whole dataset to generate the small input as well
  install_deps_source_setup $1

  mkdir -p "$outputs_dir"

  scripts_inputs=(
      "nfa-regex;10M.txt" # 100M
      "sort;100M.txt"
      "top-n;100M.txt"
      "wf;100M.txt"
      "spell;100M.txt"
      "diff;100M.txt"
      "bi-grams;100M.txt"
      "set-diff;100M.txt"
      "sort-sort;100M.txt"
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
    # source the required variables from setup.sh
    source_var ${1:-"notsmall"} $input
    printf -v pad %30s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:30}

    seq_outputs_file="${outputs_dir}/${script}.${seq_outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$seq_outputs_file"; } 2>&1) | tee -a "$seq_times_file"
  done

  cd ..
}

oneliners_pash(){
  times_file="par.res"
  outputs_suffix="par.out"
  time_suffix="par.time"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs"
  if [ -e "oneliners/$times_file" ]; then
    echo "skipping oneliners/$times_file"
    return 0
  fi
  
  cd oneliners/

  install_deps_source_setup $1

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  scripts_inputs=(
      "nfa-regex;10M.txt" # 100M.txt"
      "sort;100M.txt" # 3G.txt
      "top-n;100M.txt"
      "wf;100M.txt"
      "spell;100M.txt"
      "diff;100M.txt"
      "bi-grams;100M.txt"
      "set-diff;100M.txt"
      "sort-sort;100M.txt"
      "shortest-scripts;all_cmdsx100.txt"
  )

  touch "$times_file"
  echo executing one-liners with pash $(date) | tee -a "$times_file"
  echo '' >> "$times_file"

  width="2"

  for script_input in ${scripts_inputs[@]}
  do
    IFS=";" read -r -a script_input_parsed <<< "${script_input}"
    script="${script_input_parsed[0]}"
    input="${script_input_parsed[1]}"
    source_var $1 $input
    printf -v pad %30s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:30}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    pash_log="${pash_logs_dir}/${script}.pash.log"
    single_time_file="${outputs_dir}/${script}.${time_suffix}"

    echo -n "${padded_script}" | tee -a "$times_file"
    { time "$PASH_TOP/pa.sh" -w "${width}" $PASH_FLAGS --bash --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done

  cd ..
}

# oneliners --full
# oneliners_pash --full

# compare_outputs "oneliners/outputs"

nlp(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "nlp/${times_file}" ]; then
    echo "skipping nlp/${times_file}"
    return 0
  fi

  cd nlp/

  install_deps_source_setup $1

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
  echo executing Unix-for-nlp $(date) | tee -a "$times_file"
  echo '' >> "$times_file"

  source_var $1

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
  cd ..
}

nlp_pash(){
  times_file="par.res"
  outputs_suffix="par.out"
  time_suffix="par.time"
  outputs_dir="outputs"
  pash_logs_dir="pash_logs"
  width=16
  if [ -e "nlp/${times_file}" ]; then
    echo "skipping nlp/${times_file}"
    return 0
  fi

  cd nlp/

  install_deps_source_setup $1

  mkdir -p "$outputs_dir"
  mkdir -p "$pash_logs_dir"

  names_scripts=(
    "1syllable_words;6_4"
    "2syllable_words;6_5"
    "4letter_words;6_2"
    "bigrams_appear_twice;8.2_2"
    "bigrams;4_3"
    "compare_exodus_genesis;8.3_3"
    "count_consonant_seq;7_2"
    "count_morphs;7_1"
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
  echo executing Unix-for-nlp with pash $(date) | tee -a "$times_file"
  echo '' >> "$times_file"

  source_var $1

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
    { time "$PASH_TOP/pa.sh" -w "${width}" $PASH_FLAGS --log_file "${pash_log}" ${script}.sh > "$outputs_file"; } 2> "${single_time_file}"
    cat "${single_time_file}" | tee -a "$times_file"
  done
  cd ..
}

nlp --full
nlp_pash --full

compare_outputs "nlp/outputs"