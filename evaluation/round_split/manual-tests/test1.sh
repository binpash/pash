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
testFile=../../evaluation/scripts/input/10M.txt
batchSize=100000
mkfifo $file1
mkfifo $file2
mkfifo $file3
mkfifo $file4
# mkfifo $file5

mkfifo $file7
mkfifo $file8
mkfifo $file9

# mkfifo $file6
# cat $testFile > $file9 &
# ../auto-split.sh $file9 $file1 $file2 &
# grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file3 > $file4 &
# grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file2 > $file6 &
# ../eager.sh $file1 $file3 temp &
# ../eager.sh $file6 $file7 temp2 &
# cat $file4 $file7 > $file5

../r_split $testFile $batchSize $file1 $file2 $file7 &

../r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file1 > $file3 &
../r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file2 > $file4 &
../r_wrap grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' < $file7 > $file8 &

../r_merge $file3 $file4 $file8> $file5

# cat $testFile | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > $file6
# if cmp -s "$file6" "$file5"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file6" "$file5"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file6" "$file5"
# fi

# rm -rf *out