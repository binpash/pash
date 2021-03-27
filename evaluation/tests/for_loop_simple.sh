for _ in $times; do
    cat $IN | tr A-Z a-z | sort
done
