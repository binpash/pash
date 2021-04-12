[ $(uname) = 'Darwin' ] && a=/usr/share/dict/web2 || a=/usr/share/dict/words

if [ -f $a ]; then
  cat $a $a $a $a $a $a $a $a | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' | wc -l
else
  echo "Dictionary file $a not found.."
fi

