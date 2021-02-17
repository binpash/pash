INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/log_parsing/201011271400.dump}
tcpdump -nn -r ${INPUT} -A 'port 53'| sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' > /dev/null                    
tcpdump -nn -r ${INPUT} -A 'port 53'| sort |uniq |grep -Ev '(com|net|org|gov|mil|arpa)' > /dev/null
#tcpdump -nn -r ${INPUT} -A 'port 53'| grep -i referer
#tcpdump -nn -r ${INPUT} -A 'port 53'| uniq |grep -i referer
#tcpdump -nn -r ${INPUT} -A 'port 53'| grep -Ev '(GET|HEAD)'
#tcpdump -nn -r ${INPUT} -A 'port 53'| uniq |grep -Ev '(GET|HEAD)'
