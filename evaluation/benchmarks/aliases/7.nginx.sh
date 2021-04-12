#!/bin/bash
# fetch hit count for each ip ?
set -e
IN=${IN:-${PASH_TOP}/evaluation/benchmarks/aliases/input}
# original command  tail -10000 /var/log/nginx/access.log | cut -d "" "" -f1 | sort | uniq -c | sort -n | tail -n 30 | sort -nrk 1 | awk
cat ${IN}/access.log | cut -d ' ' -f1 | sort | uniq -c | sort -n | tail -n 30 | sort -nrk 1  
