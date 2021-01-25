file1=1.out
file2=2.out
file3=3.out
file4=4.out
file5=5.out
testFile=test1.txt
batchSize=100

../r_split $testFile $batchSize $file1 $file2

../r_wrap cat < $file1 > $file3
../r_wrap cat < $file2 > $file4 

../r_merge $file3 $file4 > $file5

if cmp -s "$testFile" "$file5"; then
    printf 'The file "%s" is the same as "%s"\n' "$testFile" "$file5"
else
    printf 'The file "%s" is different from "%s"\n' "$testFile" "$file5"
fi

# rm -rf *.out