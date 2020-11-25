#!/bin/bash
# Issue: Deadlocks as s2 is read from more than one parallel copy
# mkfifo s1 s2
# TEMP_C1="/tmp/{/}.out1"
# TEMP_C2="/tmp/{/}.out2"
# mkfifo ${TEMP1} ${TEMP2}
# parallel "cat {} | cut -d ' ' -f 1 | sort > $TEMP_C1" ::: $IN &
# sort -m ${TEMP1} > s1 &
# parallel "cat {} | cut -d ' ' -f 1 | sort > $TEMP_C2" ::: $IN &
# sort -m ${TEMP2} > s2 &
# cat s1 | parallel -k --pipe --jobs ${JOBS} --block 100K "comm -23 - s2"
# rm ${TEMP1} ${TEMP2}
# rm s1 s2

## Correct
mkfifo s1
rm -f s2
TEMP_C1="/tmp/{/}.out1"
TEMP_C2="/tmp/{/}.out2"
mkfifo ${TEMP1} ${TEMP2}
parallel "cat {} | cut -d ' ' -f 1 | sort > $TEMP_C1" ::: $IN &
sort -m ${TEMP1} > s1 &
parallel "cat {} | cut -d ' ' -f 1 | sort > $TEMP_C2" ::: $IN &
sort -m ${TEMP2} > s2
cat s1 | parallel -k --pipe --jobs ${JOBS} --block 250M "comm -23 - s2"
rm ${TEMP1} ${TEMP2}
rm s1 s2
