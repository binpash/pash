#!/bin/bash
# compress all files in a directory
pure_func() {
    tempfile=$(mktemp)
    cat > $tempfile 
    gzip -c $tempfile
    rm -f $tempfile
}

for item in $(cat "$PASH_TOP/evaluation/benchmarks/file-enc/pcap_list" | head -n ${ENTRIES}); do
    cat $IN$item | pure_func > $OUT$item.out
done

echo 'done';
