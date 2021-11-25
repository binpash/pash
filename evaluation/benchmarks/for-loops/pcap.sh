#!/bin/bash
#tag: pcap analysis
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input/pcap_data}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/input/output/pcap-analysis}
LOGS=${OUT}/logs
mkdir -p ${LOGS}
run_tests() {
    INPUT=$1
    /usr/sbin/tcpdump -nn -r ${INPUT} -A 'port 53'| sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)'
    # extract URL
    /usr/sbin/tcpdump -nn -r ${INPUT} -s 0 -v -n -l | egrep -i "POST /|GET /|Host:"
    # extract passwords
    /usr/sbin/tcpdump -nn -r ${INPUT} -s 0 -A -n -l | egrep -i "POST /|pwd=|passwd=|password=|Host:"
}
export -f run_tests

for f in ${IN}/*; do
    echo $f
    logname=$OUT/$(basename $f).log
    run_tests $f &> $logname
done

echo 'done';
