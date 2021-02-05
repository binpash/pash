IN=../scripts/input/1M.txt
sort $IN > other-input.txt
OTHER=other-input.txt
functions1=("sort" "grep king" "tr [:lower:] [:upper:]" "grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'")

for f in "${functions1[@]}"; do
    echo $f
    PER_SCRIPT="$f:"
    echo "cat $IN | $f" > test-script.sh
    for width in $(seq 6); do
        $PASH_TOP/pa.sh -w $width --log_file log.txt -d 1 test-script.sh > out.txt 
        grep "Eager" log.txt
    done
done

functions2=("diff -B" "comm -23")
for f in "${functions2[@]}"; do
    echo $f
    PER_SCRIPT="$f:"
    echo "$f $IN $OTHER" > test-script.sh
    for width in $(seq 6); do
        $PASH_TOP/pa.sh -w $width --log_file log.txt -d 1 test-script.sh > out.txt
        grep "Eager" log.txt
    done
done

microbenchmark="comm-par-test2.sh diff.sh double_sort.sh for_loop_simple.sh fun-def.sh grep.sh micro_1000.sh minimal_grep.sh minimal_sort.sh no_in_script.sh sed-test.sh set-diff.sh shortest_scipts.sh sort.sh topn.sh wf.sh"
for f in $microbenchmark; do
  echo $f
  for width in $(seq 6); do
    $PASH_TOP/pa.sh -w $width --log_file log.txt -d 1 $f > out.txt
    grep "Eager" log.txt
  done
done