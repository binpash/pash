#!/bin/bash
# cat $IN | tr A-Z a-z | sort
TEMP_C="/tmp/{/}.out"
mkfifo ${TEMP}
parallel -k "cat {} | tr A-Z a-z | sort > $TEMP_C" ::: $IN & 
sort -m ${TEMP}
rm ${TEMP}
