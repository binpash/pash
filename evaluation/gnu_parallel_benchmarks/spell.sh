## Issue 1: Very ugly. 
## Issue 2: We have to know a valid block size to use
TEMP_C1="/tmp/{/}.out1"
mkfifo ${TEMP1}
parallel "cat {} | col -bx | tr -cs A-Za-z '\n' | tr A-Z a-z | tr -d '[:punct:]' | sort > $TEMP_C1" ::: $IN &
sort -m ${TEMP1} | parallel -k --jobs ${JOBS} --pipe --block 250M "uniq" | uniq | parallel -k --jobs ${JOBS} --pipe --block 250M "comm -23 - $dict"
rm ${TEMP1}
