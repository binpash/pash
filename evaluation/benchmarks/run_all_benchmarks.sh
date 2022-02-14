#!/bin/bash

## Determines whether the experimental pash flags will be tested. 
## By default they are not.
export EXPERIMENTAL=0
export DEBUG=0

for item in $@
do
    if [ "--experimental" == "$item" ]; then
        export EXPERIMENTAL=1
    fi

    if [ "--debug" == "$item" ]; then
        export DEBUG=1
    fi
done

## This script is necessary to ensure that sourcing happens with bash
source run.seq.sh
source run.par.sh

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

if [ "$EXPERIMENTAL" -eq 1 ]; then
  export PASH_FLAGS="--r_split --dgsh_tee --r_split_batch_size 1000000"
  # --speculation quick_abort is not maintained at the moment 
else
  export PASH_FLAGS=""
fi

## Add the debug flag
if [ "$DEBUG" -eq 1 ]; then
  export PASH_FLAGS="$PASH_FLAGS -d 1"
fi


oneliners
oneliners_pash

compare_outputs "oneliners/outputs"

unix50
unix50_pash

compare_outputs "unix50/outputs"

nlp
nlp_pash

compare_outputs "nlp/outputs"

web-index
web-index_pash

compare_outputs "web-index/outputs"

analytics-mts
analytics-mts_pash

compare_outputs "analytics-mts/outputs"
