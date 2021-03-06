#!/bin/bash

# FIXME: diff below ignores whitespace
# FIXME: pass flags to AGG

if [[ -z "$PASH_TOP" ]]; then
    echo "Must set PASH_TOP" 1>&2
    exit 1
fi

IN=$PASH_TOP/evaluation/scripts/input/1M.txt

CMD="$1"
FLG="$2"
# AGG=$(grep aggregate ../$CMD.json | cut -d ':' -f '2' | awk '{$1=$1};1' | tr '"' ' ' | sed 's/,$//')
AGG=$(cat ../$CMD.json | jq -r '.aggregate' | tr -d '$' | sed "s;PASH_TOP;$PASH_TOP;")

diff -b <(cat $IN $IN | $CMD $FLG) <(python $AGG <(cat $IN | $CMD $FLG) <(cat $IN | $CMD $FLG)) | head
if [ $? -ne 0 ]; then
    echo $CMD "$FLG...fail"
else
    echo $CMD "$FLG...pass"
fi



