#!/bin/bash

# FIXME: skip running if output file exists (using tee?)

set -e

# time: print real in seconds, to simplify parsing
TIMEFORMAT="%3R" # %3U %3S"

if [[ -z "$PASH_TOP" ]]; then
  echo "Need to provide PASH_TOP, possibly $(git rev-parse --show-toplevel)" 1>&2
  exit 1
fi

oneliners(){
  cd oneliners/
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi
  
  cd ./input/
  ./setup.sh
  cd ..

  echo '' > seq.res
  echo executing one-liners $(date) | tee -a ./seq.res
  export IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1M.txt}
  echo '' > seq.res
  echo nfa-regex.sh:        $({ time ./nfa-regex.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo sort.sh:             $({ time ./sort.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo top-n.sh:            $({ time ./top-n.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo wf.sh:               $({ time ./wf.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo spell.sh:            $({ time ./spell.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo bi-grams.sh:         $({ time ./bi-grams.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo diff.sh:             $({ time ./diff.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo set-diff.sh:         $({ time ./set-diff.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo shortest-scripts.sh: $({ time ./shortest-scripts.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  echo sort-sort.sh:        $({ time ./sort-sort.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  cd ..
}

unix50(){
  cd unix50/
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  cd input/
  ./setup.sh
  cd ..

  echo '' > seq.res
  echo executing unix50 $(date) | tee -a ./seq.res
  # FIXME this is the input prefix; do we want all to be IN 
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
  cd web-index/
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

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
  cd max-temp/
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

  echo '' > seq.res
  echo executing max temp $(date) | tee -a ./seq.res
  echo mex-temp.sh: $({ time ./max-temp.sh > /dev/null; } 2>&1) | tee -a ./seq.res
  cd ..
}

analytics-mts(){
  cd analytics-mts/
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

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
  cd poets/
  if [ -e ./seq.res ]; then
    echo "skipping $(basename $(pwd))/seq.res"
    return 0
  fi

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

dgsh() {
    cd dgsh
    if [ -e ./seq.res ]; then
        echo "skipping $(basename $(pwd))/seq.res"
        return 0
    fi
    ./setup.sh
    echo '' > seq.res
    echo executing dgsh $(date) | tee -a ./seq.res

    export VOC=/usr/share/dict/words
    export IN=$PASH_TOP/evaluation/benchmarks/dgsh/input
    export OUT=$PASH_TOP/evaluation/benchmarks/dgsh/input
    export BIN=/usr/local/bin



    echo compressionbench: $({ time ./1.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo gitstats: $({ time ./2.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo cmetrics: $({ time new_scripts/3.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo dublicatefiles: $({ time ./4.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo highlightwords: $({ time ./5.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo wordproperties: $({ time news_cripts/6.sh > /dev/null; } 2>&1) | tee -a ./seq.res 
    #echo weatherreport: $({ time ./7.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo textproperties: $({ time ./8.sh > /dev/null; } 2>&1) | tee -a ./seq.res 
    echo staticsymbols: $({ time ./9.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo hierarchymap: $({ time ./10.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    #echo plotgit: $({ time ./11.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo parallelword: $({ time ./12.sh > /dev/null; } 2>&1) | tee -a ./seq.res 
    #echo venuauthor: $({ time ./13.sh > /dev/null; } 2>&1) | tee -a ./seq.res #FIXME
    #echo 2dfourier: $({ time ./14.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo nuclear: $({ time new_scripts/15.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    #echo fft: $({ time ./16.sh > /dev/null; } 2>&1) | tee -a ./seq.res
    echo reordercol: $({ time ./17.sh > /dev/null; } 2>&1) | tee -a ./seq.res
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







posh
#poets
#aliases
#dgsh
