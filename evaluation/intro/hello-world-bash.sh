#!/usr/bin/env bash

# uses bashisms to demonstrate pash --bash support

[[ $(uname) = 'Darwin' ]] && a=/usr/share/dict/web2 || a=/usr/share/dict/words

pattern=('\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4')

if [[ -f $a ]]; then
  cat "$a"{,,,,,,,,} | grep "${pattern[@]}" | wc -l
else
  echo "Dictionary file $a not found.."
fi

