curl -s "$IN_DIR/p3a.txt" | gunzip | cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1
