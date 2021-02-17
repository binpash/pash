# count the packet number in a pcap file
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/log_parsing/2019-02-28-20-50-15-192.168.1.193.pcap}
tcpdump -nn -r ${INPUT} | wc -l
