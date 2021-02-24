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

# ../auto-split.sh $testFile $file1 $file2 &
# grep 'Bell' < $file1 > $file3 &
# grep 'Bell' < $file2 > $file4 &
# ../eager.sh $file4 $file6 temp &
# cat $file3 $file6 > $file5


../../runtime/r_split $testFile $batchSize $file1 $file2 &

# ../r_wrap grep 'Bell' < $file1 > $file3 &
# ../r_wrap grep 'Bell' < $file2 > $file4 & 

../../runtime/r_unwrap < $file1 > $file3 &
../../runtime/r_unwrap < $file2 > $file4 &

wc $file3 > $file5 &
wc $file4 > $file6 &
# paste -d '+' $file5 $file6 | head -n +3 | bc  | tr -s '\n'  ' ' | sed 's/^/   /' | sed 's/$/\ /'
# sum=$(($(cat $file5 ) + $(cat $file6)))
# echo sum

# wc $testFile
# if cmp -s "$file8" "$file7"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file8" "$file7"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file8" "$file7"
# fi
./merge-wc.sh $file5 $file6 > $file7
echo "${testFile}" | cat $file7 - 


rm -rf *out