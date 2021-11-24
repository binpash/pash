set_diff()
{
    grep -vx -f $1 -
}

dict=$PASH_TOP/evaluation/tests/input/sorted_words
IN=$PASH_TOP/evaluation/tests/input/1M.txt

cat $IN |
    # groff -t -e -mandoc -Tascii |  # remove formatting commands
    col -bx |                      # remove backspaces / linefeeds
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |                   # map upper to lower case
    tr -d '[:punct:]' |            # remove punctuation
    sort |                         # put words in alphabetical order
    uniq |                         # remove duplicate words
    set_diff $dict                 # report words not in dictionary 
