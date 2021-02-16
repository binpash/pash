#!/bin/bash
file1=1.out
file2=2.out
file3=3.out
file4=4.out
file5=5.out
file6=6.out
file7=7.out
file8=8.out

IN10=../../evaluation/scripts/input/10M.txt
# cat $IN10 | grep 'Bell' | cut -f 2
testFile=../../evaluation/scripts/input/10M.txt
batchSize=1000000

mkfifo $file1
mkfifo $file2
mkfifo $file3
mkfifo $file4
mkfifo $file5
mkfifo $file6

# ../auto-split.sh $testFile $file1 $file2 &
# grep 'Bell' < $file1 > $file3 &
# grep 'Bell' < $file2 > $file4 &
# ../eager.sh $file4 $file6 temp &
# cat $file3 $file6 > $file5


../r_split $testFile $batchSize $file1 $file2 &

# ../r_wrap grep 'Bell' < $file1 > $file3 &
# ../r_wrap grep 'Bell' < $file2 > $file4 & 

../r_unwrap < $file1 > $file3 &
../r_unwrap < $file2 > $file4 &

wc $file3 > $file5 &
wc $file4 > $file6 &
../agg/opt/wc.sh $file5 $file6 > $file7
# sum=$((cat $file5 + cat $file6))
# echo sum

# wc $testFile > $file8
# if cmp -s "$file8" "$file7"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file8" "$file7"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file8" "$file7"
# fi

rm -rf *out