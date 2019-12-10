cat $IN | tail +2 | paste $IN - | sort | uniq
