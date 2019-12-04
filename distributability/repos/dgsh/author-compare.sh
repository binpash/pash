#!/usr/bin/env dgsh
#
# SYNOPSIS Venue author compare
# DESCRIPTION
# Given the specification of two publication venues, read a compressed
# DBLP computer science bibliography from the standard input (e.g. piped
# from curl -s http://dblp.uni-trier.de/xml/dblp.xml.gz or from a locally
# cached copy) and output the number of papers published in each of the
# two venues as well as the number of authors who have published only in
# the first venue, the number who have published only in the second one,
# and authors who have published in both.  The venues are specified through
# the script's first two command-line arguments as a DBLP key prefix, e.g.
# journals/acta/, conf/icse/, journals/software/, conf/iwpc/, or conf/msr/.
# Demonstrates the use of dgsh-wrap -e to have sed(1) create two output
# streams and the use of tee to copy a pair of streams into four ones.
#
#  Copyright 2017 Diomidis Spinellis
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

# Extract and sort author names
sorted_authors()
{
  sed -n 's/<author>\([^<]*\)<\/author>/\1/p' |
  sort
}

# Escape a string to make it a valid sed(1) pattern
escape()
{
  echo "$1" | sed 's/\([/\\]\)/\\\1/g'
}

export -f sorted_authors

if [ ! "$2" -a ! "$DGSH_DOT_DRAW"] ; then
  echo "Usage: $0 key1 key2" 1>&2
  echo "Example: $0 conf/icse/ journals/software/" 1>&2
  exit 1
fi

gzip -dc |
# Output the two venue authors as two output streams
dgsh-wrap -e sed -n "
/^<.*key=\"$(escape $1)/,/<title>/ w >|
/^<.*key=\"$(escape $2)/,/<title>/ w >|" |
# 2 streams in 4 streams out: venue1, venue2, venue1, venue2
tee |
{{
  {{
    echo -n "$1 papers: "
    grep -c '^<.* mdate=.* key='
    echo -n "$2 papers: "
    grep -c '^<.* mdate=.* key='
  }}

  {{
    call sorted_authors
    call sorted_authors
  }} |
  comm |
  {{
    echo -n "Authors only in $1: "
    wc -l
    echo -n "Authors only in $2: "
    wc -l
    echo -n 'Authors common in both venues: '
    wc -l
  }}
}} |
cat
