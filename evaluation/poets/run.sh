for d in scripts/* ; do
    echo "$d"
    $PASH_TOP/pa.sh $d > /dev/null
done

