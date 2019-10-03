#!/bin/bash
# https://stackoverflow.com/questions/14927895/sieve-of-eratosthenes-unix-script
limit=$1
sieve="$(seq 2 $limit|sort)"

for n in 2 $(seq 3 2 $limit)
do
  sieve="$(comm -23 <(echo "$sieve") <(seq $(($n * $n)) $n $limit|sort))"
done

echo "$sieve"|sort -n
