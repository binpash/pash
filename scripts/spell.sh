# dict depends on the system
dict=/usr/share/dict/words

cat manpage.1               |  # other tricks
groff -t -e -mandoc -Tascii |  # remove formatting commands
# prepare filename | 
col -bx |                      # remove backspaces / linefeeds
tr A-Z a-z |                   # map upper to lower case
tr -d '[:punct:]' |            # remove punctuation
sort |                         # put words in alphabetical order
unique |                       # remove duplicate words
comm -13 $dict -               # report words not in dictionary 

# The 2nd `tr` is also an example  of writing a stage in multiple ways
# some trivially parallelizable (tr), others not so much (awk)
# https://stackoverflow.com/questions/49199136/removing-punctuation-from-txt-file-in-linux-terminal-using-tr-and-awk-comman
