cat "$IN_DIR/p4.out_2_00" "$IN_DIR/p4.out_2_01" |
    cut -c 89-92 |
    grep -v 999 |
    sort -rn |
    head -n1
