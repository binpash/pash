#!/bin/bash
# Calculate mispelled words in an input
# https://dl.acm.org/doi/10.1145/3532.315102

IN=${IN:-$PASH_TOP/evaluation/benchmarks/expert-oneliners/10G.txt}

cat $IN |
    # groff -t -e -mandoc -Tascii |  # remove formatting commands
    col -bx |                      # remove backspaces / linefeeds
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |                   # map upper to lower case
    tr -d '[:punct:]' |            # remove punctuation
    sort |                         # put words in alphabetical order
    uniq |                         # remove duplicate words
    comm -23 - $dict               # report words not in dictionary 
