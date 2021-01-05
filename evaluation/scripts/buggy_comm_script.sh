#!/bin/bash

mkfifo s1 s2 s3 s4 s5

PREV_IN=../evaluation/scripts/input/1M.txt
IN=/tmp/1M.txt

cat $PREV_IN > $IN
echo "end" >> $IN

cat $IN | grep "king" | tee s4 >s3 &
comm -23 $IN s3 > s1 &
comm -23 $IN s4 > s2 &
{ ../runtime/eager s2 s5 "/tmp/eager_intermediate_#file1" & }
cat s1 s5 > /tmp/buggy.out

comm -23 <(cat $IN $IN) <(cat $IN | grep "king") > /tmp/seq.out

rm s1 s2 s3 s4 s5

diff /tmp/buggy.out /tmp/seq.out
