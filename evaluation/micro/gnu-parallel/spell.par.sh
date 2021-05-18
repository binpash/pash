## Issue 1: Very ugly. It is akward to handle multiple parallel wires. (Developer effort)
## Issue 2: We have to know a valid block size to use (Performance might vary 10K -- 27m, 250K -- 4m, 250M -- 3m)
## Issue 3: adding -k is necessary for it to be correct (if not added it is usually correct but sometimes it isn't... Difficult to spot the bug)
BLOCK_SIZE=${BLOCK_SIZE:-250M}
TEMP_C1="/tmp/{/}.out1"
TEMP1=$(seq -w 0 $(($JOBS - 1)) | sed 's+^+/tmp/in+' | sed 's/$/.out1/' | tr '\n' ' ')
TEMP1=$(echo $TEMP1)
mkfifo $TEMP1
echo $TEMP1
ls $TEMP1
echo $IN
parallel "cat {} | col -bx | tr -cs A-Za-z '\n' | tr A-Z a-z | tr -d '[:punct:]' | sort > $TEMP_C1" ::: $IN &
ls $TEMP1
echo $TEMP1
sort -m $TEMP1 # | parallel -k --jobs ${JOBS} --pipe --block "$BLOCK_SIZE" "uniq" | uniq # | parallel -k --jobs ${JOBS} --pipe --block "$BLOCK_SIZE" "comm -23 - $dict"
rm $TEMP1
