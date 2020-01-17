rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file10"
mkfifo "#file10"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file15"
mkfifo "#file15"
cat "#file10" | sort  > "#file12" &
cat "#file12" | uniq  > "#file14" &
cat "${IN}" > "#file15" &
cat "${IN}" > "#file16" &
cat "#file15" | tr -cs A-Za-z "\n" > "#file17" &
cat "#file16" | tr -cs A-Za-z "\n" > "#file18" &
cat "#file17" | tr A-Z a-z > "#file19" &
cat "#file18" | tr A-Z a-z > "#file20" &
bigram_aux_reduce "#file21" "#file22" "#file23" "#file24" "#file25" "#file26" "#file10" "#file28" "#file29" &
bigram_aux_map "#file19" "#file21" "#file22" "#file23" &
bigram_aux_map "#file20" "#file24" "#file25" "#file26" &
cat "#file14" > /tmp/distr_output/0 &
wait
rm -f "#file10"
rm -f "#file18"
rm -f "#file19"
rm -f "#file20"
rm -f "#file14"
rm -f "#file17"
rm -f "#file16"
rm -f "#file12"
rm -f "#file15"
rm -rf "/dev/shm/dish"