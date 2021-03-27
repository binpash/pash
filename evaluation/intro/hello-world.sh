if [ $(uname) = 'Darwin' ]; then
  # On OSX
  a=/usr/share/dict/web2
else
  # On Linux 
  a=/usr/share/dict/words
fi

if [ -f $a ]; then
  cat $a $a $a $a $a $a $a $a | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' | wc -l
else
  echo "Dictionary file $a not found.."
fi

