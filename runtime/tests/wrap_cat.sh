file1=1.out
file2=2.out
file3=3.out
file4=4.out
testFile=../../evaluation/scripts/input/1G.txt
batchSize=60000

mkfifo $file1
# mkfifo $file3

../r_split $testFile $batchSize $file1 &
# ../r_wrap cat < $file1 > $file3 &
../r_merge $file1 > $file4

# if cmp -s "$testFile" "$file4"; then
#     printf 'The file "%s" is the same as "%s"\n' "$file1" "$file3"
# else
#     printf 'The file "%s" is different from "%s"\n' "$file1" "$file3"
# fi

rm -rf *.out