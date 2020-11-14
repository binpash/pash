#!/bin/bash
# cat $IN | tr A-Z a-z | sort
TEMP_C1="/tmp/{/}.out1"
mkfifo ${TEMP1}
parallel "cat {} | tr A-Z a-z | sort > $TEMP_C1" ::: $IN & 
sort -m ${TEMP1}
rm ${TEMP1}
