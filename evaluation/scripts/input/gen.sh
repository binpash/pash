curl 'http://ndr.md/corpus/dummy/1M.txt' > 1M.txt

touch 1G.txt
for (( i = 0; i < 1000; i++ )); do
  cat 1M.txt >> 1G.txt
done

# more logic to truncate file --- trivial on Linux, more difficult on OS X
