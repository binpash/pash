# count the packet number in a pcap file
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/201011271400.dump}
tcpdump -nn -r ${INPUT} | wc -l
