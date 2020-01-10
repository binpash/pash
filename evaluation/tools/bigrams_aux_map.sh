#!/bin/bash

## By making tee | tail | paste its own function, we can implement it
## as a pure command separated into a generalized map and a
## reduce. Following the ParSynt work, a generalized map also keeps
## some auxiliary variables (in our case streams) to enable parallelization.

##
## Map
##

bigram_aux_map()
{
    IN=$1
    OUT=$2
    AUX_HEAD=$3
    AUX_TAIL=$4

    mkfifo s2
    mkfifo aux1
    mkfifo aux2
    cat $IN |
        tee s2 aux1 aux2 |
        tail +2 |
        paste s2 - > $OUT &

    ## The goal of this is to write the first line of $IN in the $AUX_HEAD
    ## stream and the last line of $IN in $AUX_TAIL

    ## TODO: I am not sure if using head/tail like this works or breaks
    ## the pipes
    head -n 1 aux1 > $AUX_HEAD &
    tail -n 1 aux2 > $AUX_TAIL &

    wait

    rm s2
    rm aux1
    rm aux2
}

##
## Reduce:
##

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
    
    mkfifo intermediate

    cat $AUX_HEAD1 > $AUX_HEAD_OUT &
    cat $AUX_TAIL2 > $AUX_TAIL_OUT &
    paste $AUX_TAIL1 $AUX_HEAD2 > intermediate &
    cat $IN1 intermediate $IN2 > $OUT &

    wait

    rm intermediate
}
