mkfifo s1 s2 s3 s4 s5

export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

IN="$PASH_TOP/evaluation/tests/input/1M.txt"

sorted_in="/tmp/sorted.in"

sort $IN > $sorted_in

echo " king" | tee s4 >s3 &
grep -vx -f s3 - > s1 < $sorted_in &
grep -vx -f s4 - > s2 < $sorted_in &
## The eager is essential here or after tee to ensure non-deadlocks
{ "$PASH_TOP/runtime/eager.sh" s2  s5 "/tmp/eager_intermediate_#file1" & }
cat s1 s5 > grep-f.out

echo " king" | tee s4 >s3 &
comm -13 s3 - > s1 < $sorted_in &
comm -13 s4 - > s2 < $sorted_in &
## The eager is essential here or after tee to ensure non-deadlocks
{ "$PASH_TOP/runtime/eager.sh" s2  s5 "/tmp/eager_intermediate_#file1" & }
cat s1 s5 > comm.out

rm s1 s2 s3 s4 s5

diff grep-f.out comm.out

