#!/bin/bash

# run it by piping output to /dev/null

# scripts from https://unixgame.io/
# https://github.com/psinghbh/softsec.github.io
# input files https://github.com/psinghbh/softsec.github.io/tree/master/ctf/unixgame.io/challenges
# Which join is easier: http://www.theunixschool.com/2011/08/5-different-ways-to-join-all-lines-in.html
# 1 (default) + 3 + 1 + 1 + 6 + 1 + 1 + 3 + 5 + 9 + 3 + 2 + 1 = 37 (there are 3 missing)
# missing 8.5, 9.5, 12.1 

# 1.0: extract the last name
cat 1.txt | cut -d ' ' -f 2

# 1.1: extract names and sort
cat 1.txt | cut -d ' ' -f 2 | sort

# 1.2: extract names and sort
cat 1.txt | head -n 2 | cut -d ' ' -f 2

# 1.3: sort top first names
cat 1.txt | cut -d ' ' -f 1 | sort | uniq -c | sort -r

# 2.1: get all Unix utilities
cat 2.txt | cut -d ' ' -f 4 | tr -d ','

# 3.1: get lowercase first letter of last names (awk)
cat 3.txt | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

# 4.1: find number of rounds --- fmt -w1 essentially breaks at space; can be written as tr -cs or sth
cat 4.txt | fmt -w1 | grep '\.' | wc -l

# 4.2: find pieces captured by Belle
cat 4.txt | fmt -w1 | grep 'x' | grep '\.' | wc -l

# 4.3: find pieces captured by Belle with a pawn
cat 4.txt | fmt -w1 | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]' | wc -l

# 4.4: histogram of Belle's captures (-pawns) by each type of piece
cat 4.txt | fmt -w1 | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep '[KQRBN]' | cut -c 1-1 | sort | uniq -c | sort -nr

# 4.5: 4.4 + pawns
cat 4.txt | fmt -w1 | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq -c | sort -nr

# 4.6: piece used the most by Belle (not sure you need tail sed sed rather than # sort -r | head?)
cat 4.txt | fmt -w1 | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq | tail -n 3 | sed 2d | sed 2d

# 5.1: extract hello world
cat 5.txt | grep 'print' | cut -d '"' -f 2 | cut -c 1-12

# 6.1: order the bodies by how easy it would be to land on them in Thompson's
#      Space Travel game when playing at the highest simulation scale
cat 6.txt | awk '{print $2, $0}' | sort -nr | cut -d ' ' -f 2
# paste 6.txt 6.txt | tr '\t' ' ' | tr -s ' ' | cut -d ' ' -f 2,3 | sort -nr | cut -d ' ' -f 2
# cat 6.txt | sort -nr -k 2 | cut -d ' ' -f 1

# 7.1: identify number of AT&T unix versions
cat 7.txt | cut -f 1 | grep 'AT&T' | wc -l

# 7.2: find  most frequently occurring machine [nv: no need for fmt]
cat 7.txt | cut -f 2 | sort -n | uniq -c | sort -nr | head -n 1 | fmt -w1 | tail -n 1

# 7.3: all the decades in which a unix version was released
cat 7.txt | cut -f 4 | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/

# 8.1: count unix birth-year [nv: no fmt]
cat 8.txt | fmt -w1 | grep 1969 | wc -l

# 8.2: find Bell Labs location where Dennis Ritchie had his office
#      (nv: added last awk, could be `trim`)
cat 8.txt | grep 'Bell' | awk 'length <= 45' | cut -d ',' -f 2 | awk '{$1=$1};1'

# 8.3: find names of the four people most involved with unix
cat 8.txt | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1 | head -n 1

# 8.4: find longest words without hyphens (seems wrong)
cat 8.txt | fmt -w1 | tr -c '[a-z][A-Z]' '\n' | sort | awk 'length >= 16'

# 8.5: Find second-most-freq 8-character word(s) without hyphens
cat 8.txt > /dev/null

# 9.1: extract the word PORT
cat 9.1.txt | fmt -w1 | grep '[A-Z]' | tr '[a-z]' '\n' | grep '[A-Z]' | tr -d '\n' | cut -c 1-4

# 9.2: extract the word BELL
cat 9.2.txt | cut -c 1-1 | tr -d '\n'

# 9.3: animal that used to decorate the Unix room
cat 9.3.txt | cut -c 1-2 | tr -d '\n'

# 9.4: four corners with E centered, for an "X" configuration
cat 9.4.txt | fmt -w1 | grep '"' | sed 4d | cut -d '"' -f 2 | tr -d '\n'

# 9.5: backwards running clock, in a backwards poem
cat 9.5.txt > /dev/null

# 9.6: Follow the directions for grep
cat 9.6.txt | fmt -w1 | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'

# 9.7: Four corners
cat 9.7.txt | sed 2d | sed 2d | fmt -w1 | grep '[A-Z]' | tr -c '[A-Z]' '\n' | tr -d '\n'

# 9.8: TELE-communications
cat 9.8.txt | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 2d | sed 3d | sed 4d | tr -c '[A-Z]' '\n' | tr -d '\n'

# 9.9: 
cat 9.9.txt | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 1d | sed 2d | sed 3d | sed 5d | tr -c '[A-Z]' '\n' | tr -d '\n'

# 10.1: count Turing award recipients while working at Bell Labs
cat 10.txt | sed 1d | grep 'Bell' | cut -f 2 | wc -l

# 10.2: list Turing award recipients while working at Bell Labs
cat 10.txt | sed 1d | grep 'Bell' | cut -f 2

# 10.3: extract Ritchie's username 
cat 10.txt| grep Ritchie | cut -f2,5 | fmt -w1 | cut -c1 | tr '[A-Z]' '[a-z]' | tr -d '\n'
# cat 10.txt | grep 'Bell' | cut -f 2 | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

# 11.1: year Ritchie and Thompson receive the Hamming medal
cat 11.txt | grep 'UNIX' | cut -f 1

# 11.2: most repeated first name in the list?
cat 11.txt | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d

# 12.1: transform this list of instructions such that if the snake follows the
#       new instructions top to bottom, it ends on the location of the apple.
cat 12.txt > /dev/null
