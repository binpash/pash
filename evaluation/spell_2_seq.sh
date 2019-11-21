cat manpage.1               |  # other tricks
groff -t -e -mandoc -Tascii |  # remove formatting commands
col -bx |                      # remove backspaces / linefeeds
tr A-Z a-z |                   # map upper to lower case
tr -d '[:punct:]' |            # remove punctuation
sort |                         # put words in alphabetical order
uniq |                         # remove duplicate words
comm -13 $dict -               # report words not in dictionary 
