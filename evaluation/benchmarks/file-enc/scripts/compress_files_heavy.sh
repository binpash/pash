#!/bin/bash
# compress all files in a directory
pure_func() {
    tempfile=$(mktemp)
    cat > $tempfile 
    gzip -c $tempfile
    rm -f $tempfile
}

for item in $(cat "$PASH_TOP/evaluation/benchmarks/log-analysis/pcap_heavy_list" ); do
    for j in $(seq 1 $ENTRIES); do
        cat $IN$item | pure_func > $OUT$item.$j.stdout
    done
done

echo 'done';
