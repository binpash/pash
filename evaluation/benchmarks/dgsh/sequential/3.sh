#!/bin/bash

## Note: Needs to be run on a big git repository to make sense (maybe linux)

## Initialize the necessary temporary files
file1=$(mktemp)
file2=$(mktemp)
file3=$(mktemp)
file4=$(mktemp)

find "$@" \( -name \*.c -or -name \*.h \) -type f -print0 >"$file1"

echo -n 'FNAMELEN: '

tr \\0 \\n <"$file1" |
# Remove path
sed 's|^.*/||' |
# Maintain average
awk '{s += length($1); n++} END {
    if (n>0)
        print s / n;
    else
        print 0; }'

xargs -0 /bin/cat <"$file1" >"$file2"

sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" <"$file2" |
    cpp -P >"$file3"

# Structure definitions
echo -n 'NSTRUCT: '

egrep -c 'struct[   ]*{|struct[   ]*[a-zA-Z_][a-zA-Z0-9_]*[       ]*{' <"$file3"
#}} (match preceding openings)

# Type definitions
echo -n 'NTYPEDEF: '
grep -cw typedef <"$file3"

# Use of void
echo -n 'NVOID: '
grep -cw void <"$file3"

# Use of gets
echo -n 'NGETS: '
grep -cw gets <"$file3"

# Average identifier length
echo -n 'IDLEN: '

tr -cs 'A-Za-z0-9_' '\n' <"$file3" |
sort -u |
awk '/^[A-Za-z]/ { len += length($1); n++ } END {
    if (n>0)
        print len / n;
    else
        print 0; }'

echo -n 'CHLINESCHAR: '
wc -lc  <"$file2" |
    awk '{OFS=":"; print $1, $2}'

echo -n 'NCCHAR: '
sed 's/#/@/g' <"$file2" |
cpp -traditional -P |
wc -c |
awk '{OFMT = "%.0f"; print $1/1000}'

# Number of comments
echo -n 'NCOMMENT: '
egrep -c '/\*|//' <"$file2"

# Occurences of the word Copyright
echo -n 'NCOPYRIGHT: '
grep -ci copyright <"$file2"

# C files
find "$@" -name \*.c -type f -print0 >"$file2"

# Convert to newline separation for counting
tr \\0 \\n <"$file2" >"$file3"

# Number of C files
echo -n 'NCFILE: '
wc -l <"$file3"

# Number of directories containing C files
echo -n 'NCDIR: '
sed 's,/[^/]*$,,;s,^.*/,,' <"$file3" |
sort -u |
wc -l

# C code
xargs -0 /bin/cat <"$file2" >"$file3"

# Lines and characters
echo -n 'CLINESCHAR: '
wc -lc <"$file3" |
awk '{OFS=":"; print $1, $2}'

# C code without comments and strings
sed 's/#/@/g;s/\\[\\"'\'']/@/g;s/"[^"]*"/""/g;'"s/'[^']*'/''/g" <"$file3" |
cpp -P >"$file4"

# Number of functions
echo -n 'NFUNCTION: '
grep -c '^{' <"$file4"

# Number of gotos
echo -n 'NGOTO: '
grep -cw goto <"$file4"

# Occurrences of the register keyword
echo -n 'NREGISTER: '
grep -cw register <"$file4"

# Number of macro definitions
echo -n 'NMACRO: '
grep -c '@[   ]*define[   ][   ]*[a-zA-Z_][a-zA-Z0-9_]*(' <"$file4"
# Number of include directives
echo -n 'NINCLUDE: '
grep -c '@[   ]*include' <"$file4"

# Number of constants
echo -n 'NCONST: '
grep -ohw '[0-9][x0-9][0-9a-f]*' <"$file4" | wc -l 


# Header files
echo -n 'NHFILE: '
find "$@" -name \*.h -type f |
wc -l