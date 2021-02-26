#!/bin/bash
file1=1.out
file2=2.out
file3=3.out
file4=4.out
file5=5.out
file6=6.out
file7=7.out
file8=8.out

rm -f *.out

testFile=../../evaluation/scripts/input/10M.txt
batchSize=100000
testFile="/home/ubuntu/pash/evaluation/scripts/input/10M.txt"
if [ "$#" -gt "0" ]
 then
    testFile=$1
fi
if [ "$#" -gt "1" ]; then
    batchSize=$2
fi

mkfifo $file1
mkfifo $file2
mkfifo $file3
mkfifo $file4
mkfifo $file5
mkfifo $file6


../../runtime/r_split $testFile $batchSize $file1 $file2 &


../../runtime/r_unwrap < $file1 > $file3 &
../../runtime/r_unwrap < $file2 > $file4 &

wc $file3 > $file5 &
wc $file4 > $file6 &

./merge-wc.sh $file5 $file6


rm -rf *out