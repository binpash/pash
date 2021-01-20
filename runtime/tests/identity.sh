file1=1.out
file2=2.out
file3=3.out
testFile=test1.txt
batchSize=100

../r_split $testFile $batchSize $file1 $file2
../r_merge $file1 $file2 > $file3

if cmp -s "$testFile" "$file3"; then
    printf 'The file "%s" is the same as "%s"\n' "$testFile" "$file3"
else
    printf 'The file "%s" is different from "%s"\n' "$testFile" "$file3"
fi