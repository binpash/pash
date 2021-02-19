#!/usr/bin/env dgsh
#
# SYNOPSIS Word properties
# DESCRIPTION
# Read text from the standard input and list words
# containing a two-letter palindrome, words containing
# four consonants, and words longer than 12 characters.
#
# Demonstrates the use of dgsh-compatible paste as a gather function
#
# Example:
# curl ftp://sunsite.informatik.rwth-aachen.de/pub/mirror/ibiblio/gutenberg/1/3/139/139.txt | word-properties
#
#  Copyright 2013 Diomidis Spinellis
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
mkfifo a b c d a1 b1 c1 d1 
cat a | cat > a1 &
# List two-letter palindromes
cat b | sed 's/.*\(.\)\(.\)\2\1.*/p: \1\2-\2\1/;t
    g' | cat  > b1 &

# List four consecutive consonants
cat c | sed -E 's/.*([^aeiouyAEIOUY]{4}).*/c: \1/;t
    g' | cat > c1 &

# List length of words longer than 12 characters
cat d | awk '{if (length($1) > 12) print "l:", length($1);
    else print ""}' | cat > d1 &

cat db | tr -cs a-zA-Z \\n | sort -u | tee a b c d  > /dev/null &
paste a1 b1 c1 d1  | fgrep :
#cat d1 | fgrep :
rm a b c d a1 b1 c1 d1
