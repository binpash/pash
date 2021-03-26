file1=1.out
file2=2.out
file3=3.out
file4=4.out
testFile=../../evaluation/scripts/input/10M.txt
batchSize=70000

mkfifo $file1
mkfifo $file3

## 1. TODO: Deadlocks on merge + split (true)
## 2. Increasing batchsize deadlock
## 3. Improving wrap performance

../r_split $testFile $batchSize $file1 $file3 &
# ../r_wrap cat < $file1 > $file3 &
../r_merge $file1 $file3 > $file4

# cat $testFile > $file4

# if cmp -s "$testFile" "$file4"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file1" "$file3"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file1" "$file3"
# fi

rm -rf *.out
