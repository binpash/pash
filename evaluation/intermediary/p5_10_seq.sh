cat "$IN_DIR/p4.out_10_00" "$IN_DIR/p4.out_10_01" "$IN_DIR/p4.out_10_02" "$IN_DIR/p4.out_10_03" "$IN_DIR/p4.out_10_04" "$IN_DIR/p4.out_10_05" "$IN_DIR/p4.out_10_06" "$IN_DIR/p4.out_10_07" "$IN_DIR/p4.out_10_08" "$IN_DIR/p4.out_10_09" |
    cut -c 89-92 |
    grep -v 999 |
    sort -rn |
    head -n1
