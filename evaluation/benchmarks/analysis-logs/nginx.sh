###############################
### awk not working on pash ###
###############################
# sort by reponse codes
#pash 36 sec, bash 7 sec
INPUT=${PASH_TOP}/evaluation/scripts/input/access.log
cat ${INPUT} | cut -d "\"" -f3 | cut -d ' ' -f2 | sort | uniq -c | sort -rn  > /dev/null
# awk alternative, too slow
awk '{print $9}' ${INPUT} | sort | uniq -c | sort -rn > /dev/null
# find broken links broken links
awk '($9 ~ /404/)' ${INPUT} | awk '{print $7}' | sort | uniq -c | sort -rn > /dev/null
# for 502 (bad-gateway) we can run following command:
awk '($9 ~ /502/)' ${INPUT} | awk '{print $7}' | sort | uniq -c | sort -r > /dev/null
# Who are requesting broken links (or URLs resulting in 502)
awk -F\" '($2 ~ "/wp-admin/install.php"){print $1}' ${INPUT} | awk '{print $1}' | sort | uniq -c | sort -r > /dev/null
# 404 for php files -mostly hacking attempts
awk '($9 ~ /404/)' ${INPUT} | awk -F\" '($2 ~ "^GET .*\.php")' | awk '{print $7}' | sort | uniq -c | sort -r | head -n 20 > /dev/null
##############################
# Most requested URLs ########
awk -F\" '{print $2}' ${INPUT}  | awk '{print $2}' | sort | uniq -c | sort -r > /dev/null
# Most requested URLs containing XYZ
awk -F\" '($2 ~ "ref"){print $2}' ${INPUT} | awk '{print $2}' | sort | uniq -c | sort -r > /dev/null
