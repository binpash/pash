file1=1.out
file2=2.out
file3=3.out
file4=4.out
file5=5.out


batchSize=10000
testFile="/home/ubuntu/pash/evaluation/scripts/input/100M.txt"
if ![ $1 -eq 0 ]; then
    testFile=@1
fi
if ![ $2 -eq 0 ]; then
    testFile=@2
fi


mkfifo $file1
mkfifo $file2
mkfifo $file3
mkfifo $file4

../r_split $testFile $batchSize $file1 $file2 &

../r_wrap cat < $file1 > $file3 &
../r_wrap cat < $file2 > $file4 &

../r_merge $file3 $file4 > $file5

if cmp -s "$testFile" "$file5"; then
    printf 'The file "%s" is the same as "%s"\n' "$testFile" "$file5"
else
    printf 'The file "%s" is different from "%s"\n' "$testFile" "$file5"
fi

rm -rf *.out