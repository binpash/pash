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

testFile="$PASH_TOP/evaluation/scripts/input/100M.txt"
batchSize=100000
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

../../runtime/r_split $testFile $batchSize $file1 $file2 &

../../runtime/r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file1 > $file3 &
../../runtime/r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file2 > $file4 &
# ../r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file7 > $file8 &

../../runtime/r_merge $file3 $file4

# cat $testFile | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > $file6
# if cmp -s "$file6" "$file5"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file6" "$file5"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file6" "$file5"
# fi

rm -rf *out