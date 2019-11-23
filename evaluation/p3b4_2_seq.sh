curl -s "$IN_DIR/p3a_2_1.txt" "$IN_DIR/p3a_2_2.txt" | gunzip |   cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1
