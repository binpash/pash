#!/bin/bash
#tag: pcap analysis
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/pcap_data}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/output/pcap-analysis}
LOGS=${OUT}/logs
mkdir -p ${LOGS}
run_tests() {
    INPUT=$1
    /usr/sbin/tcpdump -nn -r ${INPUT} -A 'port 53' 2> /dev/null | sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' 2> /dev/null
    # extract URL
    /usr/sbin/tcpdump -nn -r ${INPUT} -s 0 -v -n -l 2> /dev/null | egrep -i "POST /|GET /|Host:" 2> /dev/null
    # extract passwords
    /usr/sbin/tcpdump -nn -r ${INPUT} -s 0 -A -n -l 2> /dev/null | egrep -i "POST /|pwd=|passwd=|password=|Host:" 2> /dev/null
}
export -f run_tests

pkg_count=0

for item in ${IN}/*;
do
    pkg_count=$((pkg_count + 1));
    run_tests $item > ${LOGS}/${pkg_count}.log
done

echo 'done';
