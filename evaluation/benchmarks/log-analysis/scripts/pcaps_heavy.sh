#!/bin/bash

# IN=${IN:-/dependency_untangling/pcap_data}
# OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/pcap-analysis}
pure_func() {
    tempfile=$(mktemp)
    cat > $tempfile 
    tcpdump -nn -r $tempfile -A 'port 53' 2> /dev/null | sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' 2> /dev/null
    # extract URL
    tcpdump -nn -r $tempfile -s 0 -v -n -l 2> /dev/null | egrep -i "POST /|GET /|Host:" 2> /dev/null
    # extract passwords
    tcpdump -nn -r $tempfile -s 0 -A -n -l 2> /dev/null | egrep -i "POST /|pwd=|passwd=|password=|Host:" 2> /dev/null

    rm -f $tempfile
}
export -f pure_func

for item in $(seq 1 $ENTRIES); do
    for j in $(cat "$PASH_TOP/evaluation/benchmarks/log-analysis/pcap_heavy_list");do
        cat ${IN}${j} | pure_func > ${OUT}${j}.${item}.stdout.log; 
    done
done

echo 'done';