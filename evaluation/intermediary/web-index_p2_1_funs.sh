#!/bin/bash

## FIXME: This is not exactly correct, as it misses the last and first
## lines of the separated files, but it will work as a first pass.

## TODO: If the results are not good, move the uniq -c and sort in the function
wi_bigrams_aux()
{
    s2=$(mktemp -u)
    ( mkfifo $s2 > /dev/null ) ;
    tee $s2 |
        tail +2 |
        paste $s2 - |
        sort
    rm $s2
}

wi_bigram_aux_reduce()
{
    IN1=$1
    IN2=$2
    OUT=$3

    sort -m $IN1 $IN2 > $OUT
}

wi_trigrams_aux()
{
    s2=$(mktemp -u)
    s3=$(mktemp -u)
    ( mkfifo $s2 $s3 > /dev/null )

    tee $s2 |
        tail +2 |
        paste $s2 - |
        tee $s3 |
        cut -f 1 |
        tail +3 |
        paste $s3 - |
        sort

    rm $s2
}

wi_trigram_aux_reduce()
{
    IN1=$1
    IN2=$2
    OUT=$3

    sort -m $IN1 $IN2 > $OUT
}


export -f wi_bigrams_aux
export -f wi_bigram_aux_reduce
export -f wi_trigrams_aux
export -f wi_trigram_aux_reduce
