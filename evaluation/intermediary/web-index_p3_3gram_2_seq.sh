
cat "$IN_DIR/p2.out_2_00" "$IN_DIR/p2.out_2_01" |
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
