if [[ -z "$PASH_TOP" ]]; then
    echo "Must set PASH_TOP" 1>&2
    exit 1
fi

IN1=$PASH_TOP/evaluation/intro/input/100M.txt
IN2=$PASH_TOP/evaluation/intro/input/words

CMD="$1"
FLG="$2"
AGG="$3"

cat $IN1 | $CMD $FLG > partial1
cat $IN2 | $CMD $FLG > partial2
cat $IN1 $IN2 | $CMD $FLG > reference

$AGG partial1 partial2 aggregated

diff aggregated reference > log
if [ $? -ne 0 ]; then
    cat log | head
    echo $CMD "$FLG...fail"
else
    echo $CMD "$FLG...pass"
fi

rm -f partial1 partial2 aggregated reference log