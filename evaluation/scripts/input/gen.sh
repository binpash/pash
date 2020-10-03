#!/bin/bash

curl 'http://ndr.md/corpus/dummy/1M.txt' > 1M.txt

touch 10M.txt
for (( i = 0; i < 10; i++ )); do
  cat 1M.txt >> 10M.txt
done

touch 100M.txt
for (( i = 0; i < 10; i++ )); do
  cat 10M.txt >> 100M.txt
done

touch 1G.txt
for (( i = 0; i < 10; i++ )); do
  cat 100M.txt >> 1G.txt
done

# more logic to truncate file --- trivial on Linux, more difficult on OS X
