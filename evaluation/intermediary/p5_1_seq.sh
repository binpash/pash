cat "$IN_DIR/p4.out" |
    cut -c 89-92 |
    grep -v 999 |
    sort -rn |
    head -n1

