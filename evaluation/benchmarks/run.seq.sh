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
  times_file="seq.res"
  outputs_suffix="seq.out"
  outputs_dir="outputs"
  if [ -e "max-temp/${times_file}" ]; then
    echo "skipping max-temp/${times_file}"
    return 0
  fi

  cd max-temp/

  mkdir -p "$outputs_dir"

  touch "$times_file"
  echo executing max temp $(date) | tee -a "$times_file"
  outputs_file="${outputs_dir}/max-temp.${outputs_suffix}"
  echo mex-temp.sh: $({ time ./max-temp.sh > "${outputs_file}"; } 2>&1) | tee -a "$times_file"
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
  cd aliases/
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  cd meta/
  ./setup.sh
  cd ..
  echo '' > seq.res
  echo executing aliases $(date) | tee -a ./seq.res
  cd select
  export WAV=$PASH_TOP/evaluation/benchmarks/aliases/meta/wav
  export JPG=$PASH_TOP/evaluation/benchmarks/aliases/meta/jpg
  export RTF=$PASH_TOP/evaluation/benchmarks/aliases/meta/rtf
  export GIT=$PASH_TOP/evaluation/benchmarks/aliases/meta/linux
  export IN=$PASH_TOP/evaluation/benchmarks/aliases/meta/
  export OUT=$PASH_TOP/evaluation/benchmarks/aliases/meta/out
  echo tomp3: $({ time ./1.tomp3.sh > /dev/null; } 2>&1) | tee -a ../seq.res
  echo unrtf: $({ time ./2.unrtf.sh > /dev/null; } 2>&1) | tee -a ../seq.res
  echo convertjpg: $({ time ./3.resiz.sh > /dev/null; } 2>&1) | tee -a ../seq.res
  #echo gitkernel: $({ time ./4.gitkernel.sh > /dev/null; } 2>&1) | tee -a ../seq.res  FIXME need complex grep command
  echo apachelog: $({ time ./5.apachelog.sh > /dev/null; } 2>&1) | tee -a ../seq.res
  echo msg: $({ time ./6.msg.sh > /dev/null; } 2>&1) | tee -a ../seq.res 
  echo nginx: $({ time ./7.nginx.sh > /dev/null; } 2>&1) | tee -a ../seq.res
  echo varlog: $({ time ./8.varlog.sh > /dev/null; } 2>&1) | tee -a ../seq.res
}


bio() {
  cd bio
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi
  export IN=$PASH_TOP/evaluation/benchmarks/bio/
  # takes too many files to download
  export IN_N=input_all.txt
  export OUT=$PASH_TOP/evaluation/benchmarks/bio/output
  # bio4
  ./setup.sh
  echo '' > seq.res
  echo executing bio $(date) | tee -a ./seq.res

  echo bio4.sh: $({ time ./bio4.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  # echo bio2.sh: $({ time ./bio2.sh > /dev/null; } 2>&1) | tee -a ./seq.res to check
}

# everything under this line is WIP

dgsh() {
    cd dgsh
    if [ -e ./seq.res ]; then
        echo "skipping $(basename $(pwd))/seq.res"
        return 0
    fi
    ./input/setup.sh
    echo '' > seq.res
    echo executing dgsh $(date) | tee -a ./seq.res

    export VOC=/usr/share/dict/words
    export IN=$PASH_TOP/evaluation/benchmarks/dgsh/input
    export MINI=$IN/mini.xml
    export OUT=$PASH_TOP/evaluation/benchmarks/dgsh/input
    export BIN=/usr/local/bin



    echo compressionbench: $({ time ./1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo gitstats: $({ time ./2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo cmetrics: $({ time manual/3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo dublicatefiles: $({ time ./4.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo highlightwords: $({ time ./5.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo wordproperties: $({ time ./6.sh > /dev/null; } 2>&1) | tee -a ./seq.res 
    #echo weatherreport: $({ time ./7.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo textproperties: $({ time ./8.sh > /dev/null; } 2>&1) | tee -a ./seq.res 
    echo staticsymbols: $({ time ./9.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    #echo hierarchymap: $({ time ./10.sh > /dev/null; } 2>&1) | tee -a ./seq.res     #dont know how it works
    #echo plotgit: $({ time ./11.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo parallelword: $({ time ./manual/12.sh > /dev/null; } 2>&1) | tee -a ./seq.res 
    #echo venuauthor: $({ time ./13.sh > /dev/null; } 2>&1) | tee -a ./seq.res #FIXME
    #echo 2dfourier: $({ time ./14.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    #echo nuclear: $({ time ./15.sh > /dev/null; } 2>&1) | tee -a ./seq.res Cannot the trace the error, the script runs fine
    #echo fft: $({ time ./16.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo reordercol: $({ time ./manual/17.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo dirlisting: $({ time ./18.sh > /dev/null; } 2>&1) | tee -a ./seq.res
}



posh() {
    cd posh
    if [ -e ./seq.res ]; then
        echo "skipping $(basename $(pwd))/seq.res"
        return 0
    fi
    ./setup.sh
    echo '' > seq.res
    echo executing posh $(date) | tee -a ./seq.res

    #export IN=$PASH_TOP/evaluation/benchmarks/posh/input
    export OUT=$PASH_TOP/evaluation/benchmarks/posh/output

    echo discat: $({ time ./1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo convert: $({ time ./2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo raytracing: $({ time ./3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    #echo zannotate: $({ time ./4.sh > /dev/null; } 2>&1) | tee -a ./seq.res where is zannotate binary
}







#posh
#poets
#aliases
dgsh

