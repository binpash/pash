for IN in $files; do
    cat $IN $IN | tr A-Z a-z | sort
done
