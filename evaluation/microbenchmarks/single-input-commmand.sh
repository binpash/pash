IN=../scripts/input/1M.txt
DIGITS=3
TIMEFORMAT="RadhasTimer %${DIGITS}R"

functions=("seq 1 100" "wc" "sort" "cat" "grep king" "tr [:lower:] [:upper:]" "grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'" "sed 10q" "uniq -c")
PER_SCRIPT="widths:"
for width in $(seq 6); do
  echo $width > percommand.csv
  PER_SCRIPT="$PER_SCRIPT,  $(cat percommand.csv)"
done
echo $PER_SCRIPT >> ./command-timings.csv

for f in "${functions[@]}"; do
    echo $f
    PER_SCRIPT="$f:"
    echo "cat $IN | $f" > test-script.sh
    for width in $(seq 6); do
        (time $PASH_TOP/pa.sh -w $width -d 0 test-script.sh > out.txt) 2>&1 | grep RadhasTimer | sed 's/RadhasTimer //' > percommand.csv
        PER_SCRIPT="$PER_SCRIPT,  $(cat percommand.csv)"
    done
    echo $PER_SCRIPT >> ./command-timings.csv
done