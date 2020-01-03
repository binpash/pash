cat "$IN_DIR/p3a.txt" | xargs -n1 curl -s | gunzip | cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1
