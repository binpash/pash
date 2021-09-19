set -e
( { false; }
  { echo one; } ) | cat
echo two
