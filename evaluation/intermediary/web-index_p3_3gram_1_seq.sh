
cat "$IN_DIR/p3.out" |
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
