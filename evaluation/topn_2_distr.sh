rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file16"
mkfifo "#file16"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file10"
mkfifo "#file10"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file17"
mkfifo "#file17"
cat "#file10" | uniq -c > "#file12" &
cat "#file12" | sort -rn > "#file14" &
cat "#file14" | sed ${N}q > "#file16" &
cat ${IN} > "#file17" &
cat ${IN} > "#file18" &
cat "#file17" | tr -cs A-Za-z "\n" > "#file19" &
cat "#file18" | tr -cs A-Za-z "\n" > "#file20" &
cat "#file19" | tr A-Z a-z > "#file21" &
cat "#file20" | tr A-Z a-z > "#file22" &
cat "#file21" | sort  > "#file23" &
cat "#file22" | sort  > "#file24" &
sort -m --parallel=2 "#file23" "#file24" > "#file10" &
cat "#file16" > /tmp/distr_output/0 &
wait
rm -f "#file16"
rm -f "#file18"
rm -f "#file24"
rm -f "#file19"
rm -f "#file20"
rm -f "#file10"
rm -f "#file21"
rm -f "#file14"
rm -f "#file23"
rm -f "#file22"
rm -f "#file12"
rm -f "#file17"
rm -rf "/dev/shm/dish"