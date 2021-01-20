mkfifo s1 s2 s3 s4 s5

IN=../evaluation/scripts/input/1M.txt

echo " king" | tee s4 >s3 &
grep -f s3 $IN > s1 &
grep -f s4 $IN > s2 &
## The eager is essential here or after tee to ensure non-deadlocks
{ ../runtime/eager s2  s5 "/tmp/eager_intermediate_#file1" & }
cat s1 s5

rm s1 s2 s3 s4 s5

