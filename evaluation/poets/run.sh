# executes all the scripts for poets from the folder scripts/
for d in scripts/* ; do
    echo "$d"
    $PASH_TOP/pa.sh $d > /dev/null
done
