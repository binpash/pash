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

testFile=../../evaluation/scripts/input/100M.txt
batchSize=10000000
testFile="/home/ubuntu/pash/evaluation/scripts/input/100M.txt"
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
mkfifo $file7
mkfifo $file8

../../runtime/r_split $testFile $batchSize $file1 $file2 &

../../runtime/r_unwrap < $file1 > $file3 &
../../runtime/r_unwrap < $file2 > $file4 &

../../runtime/eager.sh $file3 $file5 "/tmp/pash_eager_intermediate_#file1" &
../../runtime/eager.sh $file4 $file6 "/tmp/pash_eager_intermediate_#file2" &

sort < $file5 > $file7 &
sort < $file6 > $file8 &

sort -m $file7 $file8

# cat $testFile | sort > $file8
# if cmp -s "$file7" "$file8"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file7" "$file8"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file7" "$file8"
# fi

# rm -rf *out