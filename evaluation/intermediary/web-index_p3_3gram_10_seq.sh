
cat "$IN_DIR/p2.out_10_00" "$IN_DIR/p2.out_10_01" "$IN_DIR/p2.out_10_02" "$IN_DIR/p2.out_10_03" "$IN_DIR/p2.out_10_04" "$IN_DIR/p2.out_10_05" "$IN_DIR/p2.out_10_06" "$IN_DIR/p2.out_10_07" "$IN_DIR/p2.out_10_08" "$IN_DIR/p2.out_10_09" |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    tee shift2 |
    tail +2 |
    paste shift2 - |
    tee shift3 |
    cut -f 1 |
    tail +3 |
    paste shift3 - |
    sort |
    uniq -c |
    sort -rn
