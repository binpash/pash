#!/bin/bash
file1=1.out
file2=2.out
file3=3.out
file4=4.out
file5=5.out
file6=6.out

IN10=../../evaluation/scripts/input/10M.txt
# cat $IN10 | grep 'Bell' | cut -f 2
testFile=../../evaluation/scripts/input/1G.txt
batchSize=60000

mkfifo $file1
mkfifo $file2
mkfifo $file3
mkfifo $file4
mkfifo $file6
# mkfifo $file5

# ../r_split $testFile $batchSize $file1 $file2 &

# ../r_wrap grep 'Bell' < $file1 > $file3 &
# ../r_wrap grep 'Bell' < $file2 > $file4 &

# ../r_merge $file3 $file4 > $file5


../auto-split.sh $testFile $file1 $file2 &
grep 'Bell' < $file1 > $file3 &
grep 'Bell' < $file2 > $file4 &
../eager.sh $file4 $file6 temp &
cat $file3 $file6 > $file5

# cat $testFile | grep 'Bell' > $file6
# if cmp -s "$file6" "$file5"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file6" "$file5"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file6" "$file5"
# fi
# cut -f 2 < $file5

rm -rf *out