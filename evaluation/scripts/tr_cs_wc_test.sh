## This script is used to experiment with how to get parallelism benefits from a bunch of Unix50 pipelines

## You have to run the following before running this script.
## The output should be 439M long
## Warning: Takes a long time
## cat $PASH_TOP/evaluation/unix50/4.txt | $PASH_TOP/runtime/multiply.sh -m 1000000 | pv > $PASH_TOP/evaluation/unix50/big_4.txt

FILE="${PASH_TOP}/evaluation/unix50/big_4.txt"

# cat $FILE | tr -s ' ' '\n' | grep 'x' | grep '\.' | wc -l

cat $FILE | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l

## Possible solutions:
## 1. Make an aggregator for tr -s (This is the best solutoin)
## 2. Remove the -s since it is not actually necessary
## 3. Make an aggregator for wc (?)