## Issue 1: Very ugly. It is akward to handle multiple parallel wires. (Developer effort)
## Issue 2: We have to know a valid block size to use (Performance might vary 10K -- 27m, 250K -- 4m, 250M -- 3m)
## Issue 3: adding -k is necessary for it to be correct (if not added it is usually correct but sometimes it isn't... Difficult to spot the bug)
TEMP_C1="/tmp/{/}.out1"
mkfifo ${TEMP1}
parallel "cat {} | col -bx | tr -cs A-Za-z '\n' | tr A-Z a-z | tr -d '[:punct:]' | sort > $TEMP_C1" ::: $IN &
sort -m ${TEMP1} | parallel -k --jobs ${JOBS} --pipe --block 250M "uniq" | uniq | parallel -k --jobs ${JOBS} --pipe --block 250M "comm -23 - $dict"
rm ${TEMP1}
