cat "$IN_DIR/p3a_10_1.txt" "$IN_DIR/p3a_10_2.txt" "$IN_DIR/p3a_10_3.txt" "$IN_DIR/p3a_10_4.txt" "$IN_DIR/p3a_10_5.txt" "$IN_DIR/p3a_10_6.txt" "$IN_DIR/p3a_10_7.txt" "$IN_DIR/p3a_10_8.txt" "$IN_DIR/p3a_10_9.txt" "$IN_DIR/p3a_10_10.txt" | xargs -n1 curl -s | gunzip |   cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1
