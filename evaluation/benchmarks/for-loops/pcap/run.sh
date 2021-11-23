INPUT=pcaps
OUT=logs
run_tests() {
    INPUT=$1
    echo $INPUT
    /usr/sbin/tcpdump -nn -r ${INPUT} -A 'port 53'| sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' > /dev/null                    
    # extract URL
    /usr/sbin/tcpdump -nn -r ${INPUT} -s 0 -v -n -l | egrep -i "POST /|GET /|Host:" > /dev/null
    # extract passwords
    /usr/sbin/tcpdump -nn -r ${INPUT} -s 0 -A -n -l | egrep -i "POST /|pwd=|passwd=|password=|Host:"
}
export -f run_tests
rm -rf $OUT
mkdir -p $OUT
for f in ${INPUT}/*; do
    logname=$OUT/$(basename $f .pcapng).log
    run_tests $f &> $logname
done

echo "done"
