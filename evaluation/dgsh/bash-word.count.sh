#!/usr/bin/env dgsh
#
# SYNOPSIS Parallel word count
# DESCRIPTION
# Number of processes
N=8
# Collation order for sorting
export LC_ALL=C


# split our chunk in N files
split -dn l/$N $1
# our files are in x00 x01 x02 

for i in $(seq $((N -1))) ; 
do 
    f=x0$i
    cat $f | tr -s ' \t\n\r\f' '\n' | sort -S 512M  > ${f}.res
done 

# i don't like it at all, u def can optimize it
cat *.res | sort -S 512M  | uniq -c
