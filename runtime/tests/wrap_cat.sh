file1=1.out
file2=2.out
file3=3.out
file4=4.out
testFile=test1.txt
batchSize=100

../r_split $testFile $batchSize $file1 $file2
../r_wrap cat < $file1 > $file3

if cmp -s "$file1" "$file3"; then
    printf 'The file "%s" is the same as "%s"\n' "$file1" "$file3"
else
    printf 'The file "%s" is different from "%s"\n' "$file1" "$file3"
fi

# rm -rf *.out