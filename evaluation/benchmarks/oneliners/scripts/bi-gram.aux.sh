#!/bin/bash
# Auxiliary functions for bi-grams

bigrams_aux()
{
    s2=$(mktemp -u)
    mkfifo $s2
    tee $s2 |
        tail -n +2 |
        paste $s2 - |
        sed '$d'
    rm $s2
}

bigram_aux_map()
{
    IN=$1
    OUT=$2
    AUX_HEAD=$3
    AUX_TAIL=$4

    s2=$(mktemp -u)
    aux1=$(mktemp -u)
    aux2=$(mktemp -u)
    aux3=$(mktemp -u)
    temp=$(mktemp -u)

    mkfifo $s2
    mkfifo $aux1
    mkfifo $aux2
    mkfifo $aux3

    ## New way of doing it using an intermediate file. This is slow
    ## but doesn't deadlock
    cat $IN > $temp

    sed '$d' $temp > $aux3 &
    cat $temp | head -n 1 > $AUX_HEAD &
    cat $temp | tail -n 1 > $AUX_TAIL &
    cat $temp | tail -n +2 | paste $aux3 - > $OUT &

    # ## Old way of doing it
    # cat $IN |
    #     tee $s2 $aux1 $aux2 |
    #     tail -n +2 |
    #     paste $s2 - > $OUT &

    # ## The goal of this is to write the first line of $IN in the $AUX_HEAD
    # ## stream and the last line of $IN in $AUX_TAIL

    # cat $aux1 | ( head -n 1 > $AUX_HEAD; $PASH_TOP/evaluation/tools/drain_stream.sh ) &
    # # while IFS= read -r line
    # # do
    # #     old_line=$line
    # # done < $aux2
    # # echo "$old_line" > $AUX_TAIL
    # ( tail -n 1 $aux2 > $AUX_TAIL; $PASH_TOP/evaluation/tools/drain_stream.sh ) &

    wait

    rm $temp
    rm $s2
    rm $aux1
    rm $aux2
    rm $aux3
}

bigram_aux_reduce()
{
    IN1=$1
    AUX_HEAD1=$2
    AUX_TAIL1=$3
    IN2=$4
    AUX_HEAD2=$5
    AUX_TAIL2=$6
    OUT=$7
    AUX_HEAD_OUT=$8
    AUX_TAIL_OUT=$9

    temp=$(mktemp -u)

    mkfifo $temp

    cat $AUX_HEAD1 > $AUX_HEAD_OUT &
    cat $AUX_TAIL2 > $AUX_TAIL_OUT &
    paste $AUX_TAIL1 $AUX_HEAD2 > $temp &
    cat $IN1 $temp $IN2 > $OUT &

    wait

    rm $temp
}

export -f bigrams_aux
export -f bigram_aux_map
export -f bigram_aux_reduce
