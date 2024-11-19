#!/bin/bash
# encrypt all files in a directory 

pure_func() {
    openssl enc -aes-256-cbc -k 'key'
}
export -f pure_func

for item in $(cat "$PASH_TOP/evaluation/benchmarks/file-enc/pcap_list" | head -n ${ENTRIES}); do
    cat $IN$item | pure_func > $OUT$item.out
done

echo 'done';
