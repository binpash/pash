#!/usr/bin/env dgsh
# not working
export LC_ALL=C

sort /usr/share/dict/words | cat > sort_dict
mkfifo a  sort_input
cat a |	tr -cs A-Za-z \\n | tr A-Z a-z | sort -u | cat > sort_input &
# tee process
cat $1 | tee a  > /dev/null &
comm -23 sort_input sort_dict | cat > bad && grep --fixed-strings --file=- --ignore-case --color --word-regex --context=2 < bad $1

rm -f a sort_input sort_dict bad
