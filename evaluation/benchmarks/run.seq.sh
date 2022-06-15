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
    # source the required variables from setup.sh
    source_var $1 $input
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

  install_deps_source_setup $1

  mkdir -p "$outputs_dir"

  touch "$times_file"
  echo executing Unix50 $(date) | tee -a "$times_file"
  echo '' >> "$times_file"

  source_var $1

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

  install_deps_source_setup $1

  source_var $1

  mkdir -p "$outputs_dir"
  
  touch "$times_file"
  echo executing web index $(date) | tee -a "$times_file"
  outputs_file="${outputs_dir}/web-index.${outputs_suffix}"
  echo web-index.sh: $({ time ./web-index.sh > "${outputs_file}"; } 2>&1) | tee -a "$times_file"
  cd ..
}

max-temp(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "max-temp/${times_file}" ]; then
    echo "skipping max-temp/${times_file}"
    return 0
  fi
  cd max-temp/
  
  install_deps_source_setup
  
  source_var 
  mkdir -p "$outputs_dir"
  touch "$times_file"
  echo executing max temp $(date) | tee -a "$times_file"
  outputs_file="${outputs_dir}/temp-analytics.${outputs_suffix}"
  echo max-temp.sh: $({ time ./temp-analytics.sh > "${outputs_file}"; } 2>&1) | tee -a "$times_file"
  cd ..
}

analytics-mts(){
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "analytics-mts/${times_file}" ]; then
    echo "skipping analytics-mts/${times_file}"
    return 0
  fi

  cd analytics-mts/
  install_deps_source_setup $1
 
  mkdir -p "$outputs_dir"

  touch "$times_file"
  echo executing MTS analytics $(date) | tee -a "$times_file"
  echo '' >> "$times_file"
  ## FIXME 5.sh is not working yet
  for number in `seq 4`
  do
    script="${number}"
    
    printf -v pad %20s
    padded_script="${script}.sh:${pad}"
    padded_script=${padded_script:0:20}
    # select the respective input
    source_var $1
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$times_file"
  done
  cd ..
}

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

aliases(){
  seq_times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "aliases/${seq_times_file}" ]; then
    echo "skipping aliases/${seq_times_file}"
    return 0
  fi

  cd aliases/

  cd input/
  ./setup.sh
  ./install-deps.sh
  cd ..

  mkdir -p "$outputs_dir"

  names_scripts=(
    #"tomp3;1.tomp3"
    #"unrtf;2.unrtf"
    #"convertjpg;3.resiz"
    # "gitkernel;4.gitkernel" # needs complex grep command
    "apachelog;5.apachelog"
    "msg;6.msg"
    "nginx;7.nginx"
    "varlog;8.varlog"
  )

  touch "$seq_times_file"
  echo executing aliases $(date) | tee -a "$seq_times_file"
  echo '' >> "$seq_times_file"

  export WAV=$PASH_TOP/evaluation/benchmarks/aliases/input/wav
  export JPG=$PASH_TOP/evaluation/benchmarks/aliases/input/jpg
  export RTF=$PASH_TOP/evaluation/benchmarks/aliases/input/rtf
  export GIT=$PASH_TOP/evaluation/benchmarks/aliases/input/linux
  export IN=$PASH_TOP/evaluation/benchmarks/aliases/input/
  export OUT=$PASH_TOP/evaluation/benchmarks/aliases/input/out
  for name_script in ${names_scripts[@]}
  do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"
    printf -v pad %30s
    padded_script="${name}.sh:${pad}"
    padded_script=${padded_script:0:30}

    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"

    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$seq_times_file"
  done
  cd ..
}

bio() {
  seq_times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "aliases/${seq_times_file}" ]; then
    echo "skipping aliases/${seq_times_file}"
    return 0
  fi

  cd bio/

  cd input/
  ./setup.sh
  cd ..

  mkdir -p "$outputs_dir"

  names_scripts=(
    "bio4.sh;bio4"
  )

  touch "$seq_times_file"
  echo executing bio $(date) | tee -a "$seq_times_file"
  echo '' >> "$seq_times_file"

  export IN=$PASH_TOP/evaluation/benchmarks/bio/
  # takes too many files to download
  export IN_N=input_all.txt
  export OUT=$PASH_TOP/evaluation/benchmarks/bio/output

  for name_script in ${names_scripts[@]}
  do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"
    printf -v pad %30s
    padded_script="${name}.sh:${pad}"
    padded_script=${padded_script:0:30}
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$seq_times_file"
  done
  cd ..
}

# everything under this line is WIP

dgsh() {
    seq_times_file="seq.res"
    seq_outpus_suffix="seq.out"
    outputs_dir="outputs"
    if [ -e "dgsh/$seq_times_file" ]; then
      echo "skipping dgsh/$seq_times_file"
      return 0
    fi

    cd dgsh
      
    cd ./input/
    ./setup.sh -full
    cd ..

    mkdir -p "$outputs_dir"
  
    names_scripts=(
      "compressionbench;1"
      "gitstats;2"
      "cmetrics;3"
      "dublicatefiles;4"
      "highlightwords;5"
      # "wordproperties;6"
      # "weatherreport;7"
      "textproperties;8"
      "staticsymbols;9"
      # "hierarchymap;10"
      # "plotgit;11"
      "parallelword;12"
      # "venuauthor;13"
      # "2dfourier;14"
      # "nuclear;15"
      # "fft;16"
      "reordercol;17"
      "dirlisting;18"
    )


    touch "$seq_times_file"
    echo executing DGSH $(date) | tee -a "$seq_times_file"
    echo '' >> "$seq_times_file"

    export VOC=/usr/share/dict/words
    export IN=$PASH_TOP/evaluation/benchmarks/dgsh/input
    export FULL=$IN/dblp.xml
    export MINI=$IN/mini.xml
    export OUT=$PASH_TOP/evaluation/benchmarks/dgsh/input
    export BIN=/usr/local/bin

    for name_script in ${names_scripts[@]}
    do
      IFS=";" read -r -a name_script_parsed <<< "${name_script}"
      name="${name_script_parsed[0]}"
      script="${name_script_parsed[1]}"
      printf -v pad %30s
      padded_script="${name}.sh:${pad}"
      padded_script=${padded_script:0:30}
      outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
      echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$seq_times_file"
    done
    cd ..
}

posh() {
    seq_times_file="seq.res"
    seq_outpus_suffix="seq.out"
    outputs_dir="outputs"
    if [ -e "posh/$seq_times_file" ]; then
      echo "skipping posh/$seq_times_file"
      return 0
    fi

    cd posh
      
    cd ./input/
    ./setup.sh -full
    cd ..

    mkdir -p "$outputs_dir"
  
    names_scripts=(
      "discat;1"
      "convert;2"
      "raytracing;3"
      # "zannotate;4" where is zannotate binary
    )

    touch "$seq_times_file"
    echo executing posh $(date) | tee -a "$seq_times_file"
    echo '' >> "$seq_times_file"

    export OUT=$PASH_TOP/evaluation/benchmarks/posh/input/output
    for name_script in ${names_scripts[@]}
    do
      IFS=";" read -r -a name_script_parsed <<< "${name_script}"
      name="${name_script_parsed[0]}"
      script="${name_script_parsed[1]}"
      printf -v pad %30s
      padded_script="${name}.sh:${pad}"
      padded_script=${padded_script:0:30}
      outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
      echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$seq_times_file"
    done
    cd ..
}

dependency_untangling() {
  seq_times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "dependency_untangling/${seq_times_file}" ]; then
    echo "skipping dependency_untangling/${seq_times_file}"
    return 0
  fi

  cd dependency_untangling/

  rm -rf input/output
  install_deps_source_setup $1
  mkdir -p "$outputs_dir"

  names_scripts=(
    "MediaConv1;img_convert"
    "MediaConv2;to_mp3"
    "Program_Inference;proginf"
    "LogAnalysis1;nginx"
    "LogAnalysis2;pcap"
    "Genomics_Computation;genomics"
    "AurPkg;pacaur"
    "FileEnc1;compress_files"
    "FileEnc2;encrypt_files"
  )

  touch "$seq_times_file"
  echo executing dependency_untangling $(date) | tee -a "$seq_times_file"
  echo '' >> "$seq_times_file"
  source_var 
  for name_script in ${names_scripts[@]}
  do
    IFS=";" read -r -a name_script_parsed <<< "${name_script}"
    name="${name_script_parsed[0]}"
    script="${name_script_parsed[1]}"
    printf -v pad %30s
    padded_script="${name}.sh:${pad}"
    padded_script=${padded_script:0:30}
    outputs_file="${outputs_dir}/${script}.${outputs_suffix}"
    echo "${padded_script}" $({ time ./${script}.sh > "$outputs_file"; } 2>&1) | tee -a "$seq_times_file"
  done
  cd ..
}
