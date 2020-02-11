
cat "$IN_DIR/p3.out" |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    tee shift1 |
    tail +2 |
    paste shift1 - |
    sort |
    uniq -c |
    sort -rn
