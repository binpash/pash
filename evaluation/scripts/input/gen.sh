curl 'http://ndr.md/corpus/dummy/i1M.txt' > i1M.txt

touch i1G.txt
for (( i = 0; i < 1000; i++ )); do
  cat i1M.txt >> i1G.txt
done

# more logic to truncate file --- trivial on Linux, more difficult on OS X
