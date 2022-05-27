BATCH_SIZE=$1
VIRTUAL_DIR=$2
OUTPUT1=$3
OUTPUT2=$4

tee >(
    head -n "$BATCH_SIZE" > "${VIRTUAL_DIR}/${OUTPUT1}";
    "$PASH_TOP"/evaluation/tools/drain_stream.sh &
    cat "${VIRTUAL_DIR}/${OUTPUT1}" > "${OUTPUT1}") |
    ( tail -n $((BATCH_SIZE+1)) > "${OUTPUT2}";
      "$PASH_TOP"/evaluation/tools/drain_stream.sh)
