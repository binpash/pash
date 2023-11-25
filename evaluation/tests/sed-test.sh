cat $PASH_TOP/evaluation/tests/input/1M.txt |
    sed 's;^d;da;' |
    sed 's;^;http://ndr.md/data/noaa/;' |
    sed 's;$;/;' |
    sed 's;^\(.*\)\(20[0-9][0-9]\).gz;\2/\1\2\.gz;' |
    sed 's;^;http://ndr.md/data/noaa/;' |
    sed "s#^#$WIKI#" |
    sed s/\$/'0s'/ |
    sed 1d |
    sed 4d |
    sed "\$d"