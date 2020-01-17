#!/bin/bash

alt_bigrams_aux()
{
    s2=$(mktemp -u)
    ( mkfifo $s2 > /dev/null ) ;
    ## I am not sure if this reads the stdin of the function
    tee $s2 |
        tail +2 |
        paste $s2 - |
        sort |
        uniq
    rm $s2
}

# alt_bigram_aux_map()
# {
#     IN=$1
#     OUT=$2
#     AUX_HEAD=$3
#     AUX_TAIL=$4

#     s2=$(mktemp -u)
#     aux1=$(mktemp -u)
#     aux2=$(mktemp -u)

#     mkfifo $s2
#     mkfifo $aux1
#     mkfifo $aux2
#     cat $IN |
#         tee $s2 $aux1 $aux2 |
#         tail +2 |
#         paste $s2 - |
#         sort |
#         uniq > $OUT &

#     ## The goal of this is to write the first line of $IN in the $AUX_HEAD
#     ## stream and the last line of $IN in $AUX_TAIL

#     ## TODO: I am not sure if using head/tail like this works or breaks
#     ## the pipes
#     cat $aux1 | ( head -n 1 > $AUX_HEAD; dd of=/dev/null > /dev/null 2>&1 ) &
#     tail -n 1 $aux2 > $AUX_TAIL &

#     wait

#     rm $s2
#     rm $aux1
#     rm $aux2
# }

alt_bigram_aux_reduce()
{
    IN1=$1
    IN2=$2
    OUT=$3

    sort -m $IN1 $IN2 |
        uniq > $OUT
}

export -f alt_bigrams_aux
# export -f alt_bigram_aux_map
export -f alt_bigram_aux_reduce
