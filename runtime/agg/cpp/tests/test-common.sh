CMD="$1"
FLG="$2"
AGG="$3"

cat $IN1 $IN2 | $CMD $FLG > ./temp/reference
cat $IN1 | $CMD $FLG > ./temp/partial1
cat $IN2 | $CMD $FLG > ./temp/partial2

$AGG ./temp/partial1 ./temp/partial2 $FLG > ./temp/aggregated

diff ./temp/aggregated ./temp/reference > ./temp/log
if [ $? -ne 0 ]; then
    cat ./temp/log | head
    echo $CMD "$FLG ...FAIL"
else
    echo $CMD "$FLG ...pass"
fi

rm -f ./temp/partial1 ./temp/partial2 ./temp/aggregated ./temp/reference ./temp/log
