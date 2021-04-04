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
    echo "skipping $(basename $(pwd))/$seq_times_file"
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
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  cd unix50/

  cd input/
  ./setup.sh
  cd ..

  echo '' > seq.res
  echo executing unix50 $(date) | tee -a ./seq.res
  # FIXME this is the input prefix; do we want all to be IN 
  # FIXME IN_PRE is also exported in the separate unix50 scripts, making this export here useless.
  #       I think we would like to just do it here.
  export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
  echo 1.sh:  $({ time ./1.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 2.sh:  $({ time ./2.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 3.sh:  $({ time ./3.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 4.sh:  $({ time ./4.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 5.sh:  $({ time ./5.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 6.sh:  $({ time ./6.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 7.sh:  $({ time ./7.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 8.sh:  $({ time ./8.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 9.sh:  $({ time ./9.sh  > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 10.sh: $({ time ./10.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 11.sh: $({ time ./11.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 12.sh: $({ time ./12.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 13.sh: $({ time ./13.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 14.sh: $({ time ./14.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 15.sh: $({ time ./15.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 16.sh: $({ time ./16.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 17.sh: $({ time ./17.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 18.sh: $({ time ./18.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 19.sh: $({ time ./19.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 20.sh: $({ time ./20.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 21.sh: $({ time ./21.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 22.sh: $({ time ./22.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 23.sh: $({ time ./23.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 24.sh: $({ time ./24.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 25.sh: $({ time ./25.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 26.sh: $({ time ./26.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 27.sh: $({ time ./27.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 28.sh: $({ time ./28.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 29.sh: $({ time ./29.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 30.sh: $({ time ./30.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 31.sh: $({ time ./31.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 32.sh: $({ time ./32.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 33.sh: $({ time ./33.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 34.sh: $({ time ./34.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 35.sh: $({ time ./35.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 36.sh: $({ time ./36.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  cd ..
}

web-index(){
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  cd web-index/

  cd input/
  ./setup.sh
  cd ..

  echo '' > seq.res
  echo executing web index $(date) | tee -a ./seq.res
  export IN=$PASH_TOP/evaluation/benchmarks/web-index/input/100.txt
  export WEB_INDEX_DIR=$PASH_TOP/evaluation/benchmarks/web-index/input
  export WIKI=$PASH_TOP/evaluation/benchmarks/web-index/input/
  echo web-index.sh: $({ time ./web-index.sh > /dev/null; } 2>&1) | tee -a ./seq.res
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
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  cd poets/

  cd input/
  ./setup.sh
  cd ..

  echo '' > seq.res
  echo executing Unix-for-poets $(date) | tee -a ./seq.res
  export IN=$PASH_TOP/evaluation/benchmarks/poets/input/genesis
  echo 1syllable_words.sh:               $({ time ./6_4.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 2syllable_words.sh:               $({ time ./6_5.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo 4letter_words.sh:                 $({ time ./6_2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo bigrams_appear_twice.sh:          $({ time ./8.2_2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo bigrams.sh:                       $({ time ./4_3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo compare_exodus_genesis.sh:        $({ time ./8.3_3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo count_consonant_seq.sh:           $({ time ./7_2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
# echo count_morphs.sh:                  $({ time ./7_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo count_trigrams.sh:                $({ time ./4_3b.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo count_vowel_seq.sh:               $({ time ./2_2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo count_words.sh:                   $({ time ./1_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo find_anagrams.sh:                 $({ time ./8.3_2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo merge_upper.sh:                   $({ time ./2_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo sort.sh:                          $({ time ./3_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo sort_words_by_folding.sh:         $({ time ./3_2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo sort_words_by_num_of_syllables.sh:$({ time ./8_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo sort_words_by_rhyming.sh:         $({ time ./3_3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
# echo trigram_rec.sh:                   $({ time ./6_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo uppercase_by_token.sh:            $({ time ./6_1_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo uppercase_by_type.sh:             $({ time ./6_1_2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo verses_2om_3om_2instances.sh:     $({ time ./6_7.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo vowel_sequencies_gr_1K.sh:        $({ time ./8.2_1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo words_no_vowels.sh:               $({ time ./6_3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  cd ..
}


aliases(){
  echo executing aliases
  cd aliases/select
  echo tomp3: $({ time ./1.tomp3.sh > /dev/null; } 2>&1)
  echo unrtf: $({ time ./2.unrtf.sh > /dev/null; } 2>&1)

}
