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
mkfifo $file7
mkfifo $file8

$PASH_TOP/runtime/r_split -r $testFile $batchSize $file1 $file2 &

$PASH_TOP/runtime/dgsh-tee -I -f -i $file1 -o $file5 &
$PASH_TOP/runtime/dgsh-tee -I -f -i $file2 -o $file6 &

sort < $file5 > $file7 &
sort < $file6 > $file8 &

sort -m $file7 $file8

# cat $testFile | sort > $file8
# if cmp -s "$file7" "$file8"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file7" "$file8"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file7" "$file8"
# fi

rm -rf *out