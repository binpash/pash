#!/bin/bash

# scripts from https://unixgame.io/
# https://github.com/psinghbh/softsec.github.io
# input files https://github.com/psinghbh/softsec.github.io/tree/master/ctf/unixgame.io/challenges
# Which join is easier: http://www.theunixschool.com/2011/08/5-different-ways-to-join-all-lines-in.html
# 1 (default) + 3 + 1 + 1 + 6 + 1 + 1 + 3 + 5 + 9 + 3 + 2 + 1 = 37 (there are 3 missing)
# missing 8.5, 9.5, 12.1 

if [[ -z "$IN_PRE" ]]; then
  if [[ -z "$PASH_TOP" ]]; then
    echo "Need to provide PASH_TOP, possibly $(git rev-parse --show-toplevel)" 1>&2
    exit 1
  else
    export IN_PRE=$PASH_TOP/evaluation/benchmarks/unix50/input
  fi
fi

IN1=$IN_PRE/1.txt
IN2=$IN_PRE/2.txt
IN3=$IN_PRE/3.txt
IN4=$IN_PRE/4.txt
IN5=$IN_PRE/5.txt
IN6=$IN_PRE/6.txt
IN7=$IN_PRE/7.txt
IN8=$IN_PRE/8.txt
IN91=$IN_PRE/9.1.txt
IN92=$IN_PRE/9.2.txt
IN93=$IN_PRE/9.3.txt
IN94=$IN_PRE/9.4.txt
IN95=$IN_PRE/9.5.txt
IN96=$IN_PRE/9.6.txt
IN97=$IN_PRE/9.7.txt
IN98=$IN_PRE/9.8.txt
IN99=$IN_PRE/9.9.txt
IN10=$IN_PRE/10.txt
IN11=$IN_PRE/11.txt
IN12=$IN_PRE/12.txt

# 1.0: extract the last name
cat $IN1 | cut -d ' ' -f 2

# 1.1: extract names and sort
cat $IN1 | cut -d ' ' -f 2 | sort

# 1.2: extract names and sort
cat $IN1 | head -n 2 | cut -d ' ' -f 2

# 1.3: sort top first names
cat $IN1 | cut -d ' ' -f 1 | sort | uniq -c | sort -r

# 2.1: get all Unix utilities
cat $IN2 | cut -d ' ' -f 4 | tr -d ','

# 3.1: get lowercase first letter of last names (awk)
cat $IN3 | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

# 4.1: find number of rounds
cat $IN4 | tr ' ' '\n' | grep '\.' | wc -l

# 4.2: find pieces captured by Belle
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l

# 4.3: find pieces captured by Belle with a pawn
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]' | wc -l

# 4.4: histogram of Belle's captures (-pawns) by each type of piece
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep '[KQRBN]' | cut -c 1-1 | sort | uniq -c | sort -nr

# 4.5: 4.4 + pawns
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq -c | sort -nr

# 4.6: piece used the most by Belle
cat $IN4 | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort -r | uniq | head -n 3 | tail -n 1

# 5.1: extract hello world
cat $IN5 | grep 'print' | cut -d "\"" -f 2 | cut -c 1-12

# 6.1: order the bodies by how easy it would be to land on them in Thompson's Space Travel game when playing at the highest simulation scale
cat $IN6 | awk "{print \$2, \$0}" | sort -nr | cut -d ' ' -f 2

# 7.1: identify number of AT&T unix versions
cat $IN7 | cut -f 1 | grep 'AT&T' | wc -l

# 7.2: find  most frequently occurring machine
cat $IN7 | cut -f 2 | sort -n | uniq -c | sort -nr | head -n 1 | tr -s ' ' '\n' | tail -n 1

# 7.3: all the decades in which a unix version was released
cat $IN7 | cut -f 4 | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/

# 8.1: count unix birth-year
cat $IN8 | tr ' ' '\n' | grep 1969 | wc -l

# 8.2: find Bell Labs location where Dennis Ritchie had his office
cat $IN8 | grep 'Bell' | awk 'length <= 45' | cut -d ',' -f 2 | awk "{\$1=\$1};1"

# 8.3: find names of the four people most involved with unix
cat $IN8 | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1 | head -n 1

# 8.4: find longest words without hyphens
cat $IN8 | tr -c "[a-z][A-Z]" '\n' | sort | awk "length >= 16"

# # 8.5: Find second-most-freq 8-character word(s) without hyphens
# cat $IN8 > /dev/null

# 9.1: extract the word PORT
cat $IN91 | tr ' ' '\n' | grep '[A-Z]' | tr '[a-z]' '\n' | grep '[A-Z]' | tr -d '\n' | cut -c 1-4

# 9.2: extract the word BELL
cat $IN92 | cut -c 1-1 | tr -d '\n'

# 9.3: animal that used to decorate the Unix room
cat $IN93 | cut -c 1-2 | tr -d '\n'

# 9.4: four corners with E centered, for an "X" configuration
cat $IN94 | tr ' ' '\n' | grep "\"" | sed 4d | cut -d "\"" -f 2 | tr -d '\n'

# # 9.5: backwards running clock, in a backwards poem
# cat $IN95 > /dev/null

# 9.6: Follow the directions for grep
cat $IN96 | tr ' ' '\n' | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'

# 9.7: Four corners
cat $IN97 | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'

# 9.8: TELE-communications
cat $IN98 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 2d | sed 3d | sed 4d | tr -c '[A-Z]' '\n' | tr -d '\n'

# 9.9:
cat $IN99 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 1d | sed 2d | sed 3d | sed 5d | tr -c '[A-Z]' '\n' | tr -d '\n'

# 10.1: count Turing award recipients while working at Bell Labs
cat $IN10 | sed 1d | grep 'Bell' | cut -f 2 | wc -l

# 10.2: list Turing award recipients while working at Bell Labs
cat $IN10 | sed 1d | grep 'Bell' | cut -f 2

# 10.3: extract Ritchie's username
cat $IN10 | grep 'Bell' | cut -f 2 | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

# 11.1: year Ritchie and Thompson receive the Hamming medal
cat $IN11 | grep 'UNIX' | cut -f 1

# 11.2: most repeated first name in the list?
cat $IN11 | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d


# # 12.1: transform this list of instructions such that if the snake follows the
# #       new instructions top to bottom, it ends on the location of the apple.
# cat $IN12 > /dev/null
