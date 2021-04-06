#!/bin/bash

# https://unix.stackexchange.com/questions/193441/how-can-i-implement-a-circular-flow-of-data-among-interconnected-commands

echo 1 >file

rm s1 
mkfifo s1
tail -f file |
  sed -u 's/^/1 + /'  |
  tee -a s1 > /dev/null &

cat s1 |
  xargs -0 -n 1 -d '\n' expr |
  tee -a file
