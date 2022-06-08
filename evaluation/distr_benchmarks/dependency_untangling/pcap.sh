#!/bin/bash
#tag: pcap analysis
IN=${IN:-/dependency_untangling/pcap_data}
OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/pcap-analysis}
mkdir -p $OUT

pure_func() {
    tempfile=$(mktemp)

    tee $tempfile | tcpdump -nn -r '-' -A 'port 53' 2> /dev/null | sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' 2> /dev/null
    # extract URL
    tcpdump -nn -r $tempfile -s 0 -v -n -l 2> /dev/null | egrep -i "POST /|GET /|Host:" 2> /dev/null
    # extract passwords
    tcpdump -nn -r $tempfile -s 0 -A -n -l 2> /dev/null | egrep -i "POST /|pwd=|passwd=|password=|Host:" 2> /dev/null

    rm -f $tempfile
}
export -f pure_func

for item in $(hdfs dfs -ls -C ${IN});
do
    logname=$OUT/$(basename $item).log;
    hdfs dfs -cat $item | pure_func > $logname
done

echo 'done';
