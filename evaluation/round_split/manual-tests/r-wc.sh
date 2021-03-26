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

testFile=$PASH_TOP/evaluation/scripts/input/1G.txt
batchSize=10000000
testFile="$PASH_TOP/evaluation/scripts/input/1G.txt"
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


$PASH_TOP/runtime/r_split -r $testFile $batchSize $file3 $file4 &


# $PASH_TOP/runtime/r_unwrap < $file1 > $file3 &
# $PASH_TOP/runtime/r_unwrap < $file2 > $file4 &

wc $file3 > $file5 &
wc $file4 > $file6 &

./merge-wc.sh $file5 $file6


rm -rf *out