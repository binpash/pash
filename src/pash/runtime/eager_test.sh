#!/usr/bin/env bash

mkfifo s1 s2

# IN=test_in.txt
IN=../scripts/input/1G.txt

cat "$IN" > s1 &
cat s2 > test_out.txt &
./eager s1 s2 intermediate &

wait

rm s1 s2

# diff -s $IN test_out.txt
