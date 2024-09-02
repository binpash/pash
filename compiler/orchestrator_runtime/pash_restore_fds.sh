#!/bin/bash

fds_file="${1?fd file not given}"
fdout_dir="${2?fd out dir not given}"

for outfile in $(ls -1 ${fdout_dir}); do
    cmd="cat ${fdout_dir}/${outfile} >&${outfile}"
    eval $cmd
done

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
${RUNTIME_LIBRARY_DIR}/fd_util -k -f ${fds_file}
unset fds_file
