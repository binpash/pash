INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/201011271400.dump}
INPUT2=${INPUT2:-$PASH_TOP/evaluation/scripts/input/2018-07-20-17-31-20-192.168.100.108.pcap}
tcpdump -nn -r ${INPUT} -A 'port 53'| sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' > /dev/null                    
tcpdump -nn -r ${INPUT} -A 'port 53'| sort |uniq |grep -Ev '(com|net|org|gov|mil|arpa)' > /dev/null
# without the pipes, bash takes 11 sec, with pipes, it takes 12 sec, same performance
# with pash
time tcpdump -nn -r ${INPUT2} -A -c 1000000 > /dev/null
time tcpdump -nn -r ${INPUT2} -A -c 1000000 | sort |uniq |grep -Ev '(com|net|org|gov|mil|arpa)' > /dev/null
