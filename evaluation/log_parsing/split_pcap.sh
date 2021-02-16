#!/bin/sh
# To process large pcap file, usually it is better to split it into small chunks first, 
# then process every chunk in parallel.
INPUT=${INPUT:-inputs/sample.pcap}
OUTPUT=${OUTPUT:-outputs/out.pcap}
split_size=1000
output_index=1
loop_count=10
exit_flag=0

command() {
	echo "$1" "$2" 
}

tcpdump -r ${INPUT} -w ${OUTPUT} -C ${split_size}

command ${OUTPUT}

while :
do
	loop_index=0
	while test ${loop_index} -lt ${loop_count}
	do
		if test -e ${OUTPUT}${output_index}
		then
			command ${OUTPUT} ${output_index} 
			output_index=$((output_index + 1))
			loop_index=$((loop_index + 1))
		else
			exit_flag=1
			break
		fi
	done
	wait

	if test ${exit_flag} -eq 1
	then
		exit 0
	fi
done
