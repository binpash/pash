# count the packet number in a pcap file
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input.pcap}
tcpdump -nn -r ${INPUT} | wc -l
