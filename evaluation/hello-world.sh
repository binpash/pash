if [ $(uname) = 'Darwin' ]; then
  # On OSX
  a=/usr/share/dict/web2
else
  # On Linux 
  a=/usr/share/dict/words
fi

if [ -f $a ]; then
  cat $a $a $a $a $a $a $a $a | grep -m5 -xiE '([a-z]*([a-z])\2[a-z]*){2}hello world'
else
  echo "Dictionary file $a not found.."
fi

