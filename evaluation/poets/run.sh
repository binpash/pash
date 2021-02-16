for d in scripts/* ; do
    echo "$d"
    /pash/pa.sh $d > /dev/null
done

