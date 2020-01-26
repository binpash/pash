
cat "$IN_DIR/p2.out_2_00" "$IN_DIR/p2.out_2_01" |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    tee shift1 |
    tail +2 |
    paste shift1 - |
    sort |
    uniq -c |
    sort -rn
