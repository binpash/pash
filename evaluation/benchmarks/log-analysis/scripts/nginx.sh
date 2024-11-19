#!/bin/bash
# tag: nginx logs
# IN=${IN:-/dependency_untangling/log_data}
# OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/nginx-logs}

pure_func() {
    tempfile=$(mktemp)

    tee $tempfile | cut -d "\"" -f3 | cut -d ' ' -f2 | sort | uniq -c | sort -rn   
    # awk alternative, too slow
    awk '{print $9}' $tempfile | sort | uniq -c | sort -rn  
    # find broken links broken links
    awk '($9 ~ /404/)' $tempfile | awk '{print $7}' | sort | uniq -c | sort -rn  
    # for 502 (bad-gateway) we can run following command:
    awk '($9 ~ /502/)' $tempfile | awk '{print $7}' | sort | uniq -c | sort -r  
    # Who are requesting broken links (or URLs resulting in 502)
    awk -F\" '($2 ~ "/wp-admin/install.php"){print $1}' $tempfile | awk '{print $1}' | sort | uniq -c | sort -r  
    # 404 for php files -mostly hacking attempts
    awk '($9 ~ /404/)' $tempfile | awk -F\" '($2 ~ "^GET .*.php")' | awk '{print $7}' | sort | uniq -c | sort -r | head -n 20  
    ##############################
    # Most requested URLs ########
    awk -F\" '{print $2}' $tempfile  | awk '{print $2}' | sort | uniq -c | sort -r  
    # Most requested URLs containing XYZ
    awk -F\" '($2 ~ "ref"){print $2}' $tempfile | awk '{print $2}' | sort | uniq -c | sort -r

    rm $tempfile
}
export -f pure_func

for log in $(cat "$PASH_TOP/evaluation/benchmarks/log-analysis/log_list" | head -n ${ENTRIES}); do
    cat $IN$log | pure_func > $OUT$log.out
done

echo 'done';
