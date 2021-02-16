#extract the real part
process_pkt()
(
    input=$1
    filter=$2
    show_ip_only=$3
    show_port=$4
    sort_res=$5
    uniq=$6
    extra=$7
    name=$8
    # the command we will execute, extracts both ip:port
    command="tcpdump -nn -r $input $filter"
    # extract the IP from the packet
    if [ "$show_ip_only" == "1" ]; then
        command+="| cut -f 3 -d \" \" "
        # remove port and extract only the ip
        if [ "$show_port" == "0" ]; then
            command+="| cut -f 1-4 -d \".\" "
        fi
    fi
    # sort the results
    if [ "$show_res" == "1" ]; then
        command+="| sort "
    fi
    # get unique entries
    if [ "$uniq" == "1" ]; then
        command+="| uniq "
    fi
    command+=$7
    # modify the time output to print our command
    TIMEFORMAT="$8 %Rs"
    # execute the command
    res=$(eval "time ($command>/dev/null) 2>.tmp; cat .tmp | tail -n 1")
    echo $res
)


# remove function
# the filters to apply
array="port-53 dst-port-80 tcp udp ip"
echo "$array" | tr ' ' '\n' | while read j; do
    i=$(echo \'$j\' | tr '-' ' ')
    # extract IP:port, don't sort don't, don't skip dublicates
    #process_pkt $PASH_TOP/evaluation/scripts/input.pcap "$i" 1 0 0 0 " " ip_port_$j
    ## extract IP, don't sort don't, don't skip dublicates
    #process_pkt $PASH_TOP/evaluation/scripts/input.pcap "$i" 1 1 0 0 " " ip_$j
    ## extract packet, sort, don't skip dublicates 
    #process_pkt $PASH_TOP/evaluation/scripts/input.pcap "$i" 0 0 1 0 " " sort_$j
    ## extract packet, sort, skip dublicates
    #process_pkt $PASH_TOP/evaluation/scripts/input.pcap "$i" 0 0 1 1 " " sort_uniq_$j
    #######################################################################
    ## Apply custom rules to the output
    ## run on DNS, extract packet, sort, don't skip dublicates
    process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 0 "|grep  -Ev '(com|net|org|gov|mil|arpa)' " dns_extract_sort_$j
    ## run on DNS, extract packet, sort, skip dublicates
    process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 1 "|grep  -Ev '(com|net|org|gov|mil|arpa)' " dns_extract_sort_uniq_$j
    ## extract http content
    ## run on HTTP, extract packet, sort, don't skip dublicates, get referer
    process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 0 "|grep  -i referer" http_extract_sort_referer_$j
    ## run on HTTP, extract packet, sort, skip dublicates, get referer
    process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 1 "|grep  -i referer" http_extract_sort_uniq_referer_$j
    ## run on HTTP, extract packet, sort, don't skip dublicates, GET HEAD
    process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 0 "|grep  -Ev '(GET|HEAD)'" http_extract_sort_get_head_$j
    ## run on HTTP, extract packet, sort, skip dublicates, GET HEAD
    process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 1 "|grep  -Ev '(GET|HEAD)'" http_extract_sort_uniq_get_head_$j
    ## run on HTTP, get user agent and sort unique
    #process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 0 "|grep  -Ei 'user-agent'" http_extract_sort_user-agent_$j
    ## run on HTTP, get user agent and sort unique
    #process_pkt $PASH_TOP/evaluation/scripts/input.pcap "-A $i" 0 0 1 1 "|grep  -Ei 'user-agent'" http_extract_sort_uniq_user-agent_$j
done
