#!/bin/bash

fds_file="${1?fd file not given}"

while read -r fd mode offset fpath; do
    if [ "$fd" -lt "3" ]; then
	continue
    fi
    if [ "$mode" == "r" ]; then
	cmd="exec ${fd}< $fpath"
	eval $cmd
    fi
    if [ "$mode" == "w" ]; then
	cmd="exec ${fd}>> $fpath"
	eval $cmd
    fi
    if [ "$mode" == "d" ]; then
	cmd="exec ${fd}<&${fpath}"
	eval $cmd
    fi
done < ${fds_file}
${RUNTIME_DIR}/fd_util -k -f ${fds_file}
