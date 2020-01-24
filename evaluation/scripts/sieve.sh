#!/bin/bash


# Doug McIlroy's implementation of Sieve of Eratosthenes

# A  combination of:
# https://swtch.com/~rsc/thread/
# https://stackoverflow.com/questions/14927895/sieve-of-eratosthenes-unix-script

OUT=./output/out.txt

limit=10000
sieve="$(seq 2 $limit | sort)"

for n in 2 $(seq 3 2 $limit)
do
  sieve="$(comm -23 <(echo "$sieve") <(seq $(($n * $n)) $n $limit|sort))"
done

echo "$sieve" | sort -n > $OUT
