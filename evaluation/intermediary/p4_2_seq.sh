cat "$IN_DIR/p3.out_2_00" "$IN_DIR/p3.out_2_01" | xargs -n1 curl -s | gunzip
