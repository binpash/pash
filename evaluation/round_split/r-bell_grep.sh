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

testFile="/home/tamlu/pash/evaluation/scripts/input/1G.txt"
batchSize=10000000
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

# mkfifo $file7
# mkfifo $file8
# mkfifo $file9


../../runtime/r_split $testFile $batchSize $file1 $file2 &

../../runtime/r_wrap tr A-Z a-z < $file1 > $file5 &
../../runtime/r_wrap tr A-Z a-z < $file2 > $file6 &

# ../../runtime/r_wrap grep 'Bell' < $file1 > $file5 &
# ../../runtime/r_wrap grep 'Bell' < $file2 > $file6 &
# ../r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file7 > $file8 &

../../runtime/r_merge $file5 $file6

# cat $testFile | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > $file6
# if cmp -s "$file6" "$file5"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file6" "$file5"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file6" "$file5"
# fi

rm -rf *out