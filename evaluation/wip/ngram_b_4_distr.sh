rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file19"
mkfifo "#file19"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file24"
mkfifo "#file24"
cat "#file12" | uniq  > "#file14" &
cat ${IN} > "#file15" &
cat ${IN} > "#file16" &
cat ${IN} > "#file17" &
cat ${IN} > "#file18" &
cat "#file15" | tail +2 > "#file19" &
cat "#file16" | tail +2 > "#file20" &
cat "#file17" | tail +2 > "#file21" &
cat "#file18" | tail +2 > "#file22" &
cat "#file19" | paste ${IN2} - > "#file23" &
cat "#file20" | paste ${IN2} - > "#file24" &
cat "#file21" | paste ${IN2} - > "#file25" &
cat "#file22" | paste ${IN2} - > "#file26" &
cat "#file23" | sort  > "#file27" &
cat "#file24" | sort  > "#file28" &
cat "#file25" | sort  > "#file29" &
cat "#file26" | sort  > "#file30" &
sort -m --parallel=4 "#file27" "#file28" "#file29" "#file30" > "#file12" &
cat "#file14" > /tmp/distr_output/0 &
wait
rm -f "#file19"
rm -f "#file16"
rm -f "#file30"
rm -f "#file22"
rm -f "#file18"
rm -f "#file21"
rm -f "#file29"
rm -f "#file20"
rm -f "#file17"
rm -f "#file15"
rm -f "#file28"
rm -f "#file12"
rm -f "#file26"
rm -f "#file25"
rm -f "#file14"
rm -f "#file23"
rm -f "#file27"
rm -f "#file24"
rm -rf "/dev/shm/dish"