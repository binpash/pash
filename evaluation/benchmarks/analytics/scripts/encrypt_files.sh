#!/bin/bash
# tag: file encryption

pure_func() {
    openssl enc -aes-256-cbc -pbkdf2 -iter 20000 -k 'key' -S 1234567890abcdef
}
export -f pure_func

for item in $(seq 1 "$ENTRIES"); do
    for pcap in $(cat "$PASH_TOP/evaluation/benchmarks/analytics/pcap_heavy_list"); do
        cat "${IN}${pcap}" | pure_func > "${OUT}${pcap}.${item}.enc"
    done
done

echo 'done';
