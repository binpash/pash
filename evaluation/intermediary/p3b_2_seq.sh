cat "$IN_DIR/p3a_2_1.txt" "$IN_DIR/p3a_2_2.txt" | xargs -n1 curl -s | gunzip
