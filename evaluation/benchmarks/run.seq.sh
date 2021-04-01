#!/bin/bash

set -e

# time: print real in seconds, to simplify parsing
TIMEFORMAT="%3R" # %3U %3S"

if [[ -z "$PASH_TOP" ]]; then
  echo "Need to provide PASH_TOP, possibly $(git rev-parse --show-toplevel)" 1>&2
  exit 1
fi

oneliners(){
  echo executing one-liners
  export IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1M.txt}
  cd oneliners/input
  ./setup.sh
  cd ..
  echo nfa-regex.sh:        $({ time ./nfa-regex.sh > /dev/null; } 2>&1)
  echo sort.sh:             $({ time ./sort.sh > /dev/null; } 2>&1)
  echo top-n.sh:            $({ time ./top-n.sh > /dev/null; } 2>&1)
  echo wf.sh:               $({ time ./wf.sh > /dev/null; } 2>&1)
  echo spell.sh:            $({ time ./spell.sh > /dev/null; } 2>&1)
  echo bi-grams.sh:         $({ time ./bi-grams.sh > /dev/null; } 2>&1)
  echo diff.sh:             $({ time ./diff.sh > /dev/null; } 2>&1)
  echo set-diff.sh:         $({ time ./set-diff.sh > /dev/null; } 2>&1)
  echo shortest-scripts.sh: $({ time ./shortest-scripts.sh > /dev/null; } 2>&1)
  echo sort-sort.sh:        $({ time ./sort-sort.sh > /dev/null; } 2>&1)
  cd ..
}

unix50(){
  echo executing unix50
  # FIXME this is the input prefix; do we want all to be IN 
  export IN_PRE=${IN_PRE:-$PASH_TOP/evaluation/benchmarks/unix50/input}
  cd unix50/input
  ./setup.sh
  cd ..
  echo 1.sh:  $({ time ./1.sh  > /dev/null; } 2>&1)
  echo 2.sh:  $({ time ./2.sh  > /dev/null; } 2>&1)
  echo 3.sh:  $({ time ./3.sh  > /dev/null; } 2>&1)
  echo 4.sh:  $({ time ./4.sh  > /dev/null; } 2>&1)
  echo 5.sh:  $({ time ./5.sh  > /dev/null; } 2>&1)
  echo 6.sh:  $({ time ./6.sh  > /dev/null; } 2>&1)
  echo 7.sh:  $({ time ./7.sh  > /dev/null; } 2>&1)
  echo 8.sh:  $({ time ./8.sh  > /dev/null; } 2>&1)
  echo 9.sh:  $({ time ./9.sh  > /dev/null; } 2>&1)
  echo 10.sh: $({ time ./10.sh > /dev/null; } 2>&1)
  echo 11.sh: $({ time ./11.sh > /dev/null; } 2>&1)
  echo 12.sh: $({ time ./12.sh > /dev/null; } 2>&1)
  echo 13.sh: $({ time ./13.sh > /dev/null; } 2>&1)
  echo 14.sh: $({ time ./14.sh > /dev/null; } 2>&1)
  echo 15.sh: $({ time ./15.sh > /dev/null; } 2>&1)
  echo 16.sh: $({ time ./16.sh > /dev/null; } 2>&1)
  echo 17.sh: $({ time ./17.sh > /dev/null; } 2>&1)
  echo 18.sh: $({ time ./18.sh > /dev/null; } 2>&1)
  echo 19.sh: $({ time ./19.sh > /dev/null; } 2>&1)
  echo 20.sh: $({ time ./20.sh > /dev/null; } 2>&1)
  echo 21.sh: $({ time ./21.sh > /dev/null; } 2>&1)
  echo 22.sh: $({ time ./22.sh > /dev/null; } 2>&1)
  echo 23.sh: $({ time ./23.sh > /dev/null; } 2>&1)
  echo 24.sh: $({ time ./24.sh > /dev/null; } 2>&1)
  echo 25.sh: $({ time ./25.sh > /dev/null; } 2>&1)
  echo 26.sh: $({ time ./26.sh > /dev/null; } 2>&1)
  echo 27.sh: $({ time ./27.sh > /dev/null; } 2>&1)
  echo 28.sh: $({ time ./28.sh > /dev/null; } 2>&1)
  echo 29.sh: $({ time ./29.sh > /dev/null; } 2>&1)
  echo 30.sh: $({ time ./30.sh > /dev/null; } 2>&1)
  echo 31.sh: $({ time ./31.sh > /dev/null; } 2>&1)
  echo 32.sh: $({ time ./32.sh > /dev/null; } 2>&1)
  echo 33.sh: $({ time ./33.sh > /dev/null; } 2>&1)
  echo 34.sh: $({ time ./34.sh > /dev/null; } 2>&1)
  echo 35.sh: $({ time ./35.sh > /dev/null; } 2>&1)
  echo 36.sh: $({ time ./36.sh > /dev/null; } 2>&1)
  cd ..
}

web-index(){
  echo executing web index
  cd web-index/input
  ./setup.sh
  cd ..
  export IN=$PASH_TOP/evaluation/benchmarks/web-index/input/100.txt
  export WEB_INDEX_DIR=$PASH_TOP/evaluation/benchmarks/web-index/input
  export WIKI=$PASH_TOP/evaluation/benchmarks/web-index/input/
  echo web-index.sh: $({ time ./web-index.sh > /dev/null; } 2>&1)
  cd ..
}

max-temp(){
  echo executing max temp
  cd max-temp/
  echo mex-temp.sh: $({ time ./max-temp.sh > /dev/null; } 2>&1)
  cd ..
}

analytics-mts(){
  echo executing MTS analytics
  cd analytics-mts/input
  ./setup.sh
  cd ..
  echo 1.sh: $({ time ./1.sh > /dev/null; } 2>&1)
  echo 2.sh: $({ time ./2.sh > /dev/null; } 2>&1)
  echo 3.sh: $({ time ./3.sh > /dev/null; } 2>&1)
  echo 4.sh: $({ time ./4.sh > /dev/null; } 2>&1)
#FIXME: echo 5.sh: $({ time ./5.sh > /dev/null; } 2>&1)
  cd ..
}

poets(){
  echo executing Unix-for-poets
  export IN=$PASH_TOP/evaluation/benchmarks/poets/input/pg/01frd10.txt
  cd poets/input
  ./setup.sh
  cd ..
  echo 1syllable_words.sh:               $({ time ./1syllable_words.sh > /dev/null; } 2>&1)
  echo 2syllable_words.sh:               $({ time ./2syllable_words.sh > /dev/null; } 2>&1)
  echo 4letter_words.sh:                 $({ time ./4letter_words.sh > /dev/null; } 2>&1)
  echo bigrams_appear_twice.sh:          $({ time ./bigrams_appear_twice.sh > /dev/null; } 2>&1)
  echo bigrams.sh:                       $({ time ./bigrams.sh > /dev/null; } 2>&1)
  echo compare_exodus_genesis.sh:        $({ time ./compare_exodus_genesis.sh > /dev/null; } 2>&1)
  echo count_consonant_seq.sh:           $({ time ./count_consonant_seq.sh > /dev/null; } 2>&1)
# FIXME echo count_morphs.sh:                  $({ time ./count_morphs.sh > /dev/null; } 2>&1)
  echo count_trigrams.sh:                $({ time ./count_trigrams.sh > /dev/null; } 2>&1)
  echo count_vowel_seq.sh:               $({ time ./count_vowel_seq.sh > /dev/null; } 2>&1)
  echo count_words.sh:                   $({ time ./count_words.sh > /dev/null; } 2>&1)
  echo find_anagrams.sh:                 $({ time ./find_anagrams.sh > /dev/null; } 2>&1)
  echo merge_upper.sh:                   $({ time ./merge_upper.sh > /dev/null; } 2>&1)
  echo sort.sh:                          $({ time ./sort.sh > /dev/null; } 2>&1)
  echo sort_words_by_folding.sh:         $({ time ./sort_words_by_folding.sh > /dev/null; } 2>&1)
  echo sort_words_by_num_of_syllables.sh:$({ time ./sort_words_by_num_of_syllables.sh > /dev/null; } 2>&1)
  echo sort_words_by_rhyming.sh:         $({ time ./sort_words_by_rhyming.sh > /dev/null; } 2>&1)
#FIXME echo trigram_rec.sh:                   $({ time ./trigram_rec.sh > /dev/null; } 2>&1)
  echo uppercase_by_token.sh:            $({ time ./uppercase_by_token.sh > /dev/null; } 2>&1)
  echo uppercase_by_type.sh:             $({ time ./uppercase_by_type.sh > /dev/null; } 2>&1)
  echo verses_2om_3om_2instances.sh:     $({ time ./verses_2om_3om_2instances.sh > /dev/null; } 2>&1)
  echo vowel_sequencies_gr_1K.sh:        $({ time ./vowel_sequencies_gr_1K.sh > /dev/null; } 2>&1)
  echo words_no_vowels.sh:               $({ time ./words_no_vowels.sh > /dev/null; } 2>&1)
  cd ..
}

poets
