
cat "$IN_DIR/p2.out_2_00" "$IN_DIR/p2.out_2_01" |
    sort |
    uniq -c |
    sort -rn

