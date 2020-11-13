## Issue 1: Very ugly. 
## Issue 2: We have to know a valid block size to use
TEMP_C1="/tmp/{/}.out1"
TEMP_C2="/tmp/{/}.out2"
mkfifo ${TEMP1} ${TEMP2}
parallel "cat {} | col -bx | tr -cs A-Za-z '\n' | tr A-Z a-z | tr -d '[:punct:]' | sort > $TEMP_C1" ::: $IN &
sort -m ${TEMP1} | parallel --pipe --block 250M "uniq > $TEMP_C2" &
uniq ${TEMP2} | parallel -k --pipe --block 250M "comm -23 - $dict"
rm ${TEMP1} ${TEMP2} 
