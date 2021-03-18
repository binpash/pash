#!/bin/bash
file1=1.out
file2=2.out
file3=3.out
file4=4.out
file5=5.out
file6=6.out
file7=7.out
file8=8.out
file9=9.out
rm -f *.out

testFile="$PASH_TOP/evaluation/scripts/input/10M.txt"
batchSize=1000000
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

$PASH_TOP/runtime/r_split $testFile $batchSize $file1 $file2 &

$PASH_TOP/runtime/r_wrap tr A-Z a-z < $file1 > $file3 &
$PASH_TOP/runtime/r_wrap tr A-Z a-z < $file2 > $file4 &

$PASH_TOP/runtime/r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file3 > $file5 &
$PASH_TOP/runtime/r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file4 > $file6 &

$PASH_TOP/runtime/r_merge $file5 $file6
# cat $testFile | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > t2.out
# if cmp -s t1.out t2.out; then
#     printf 'The file "%s" is the same as "%s"\n' "$file6" "$file5"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file6" "$file5"
# fi

rm -rf *out