rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file21"
mkfifo "#file21"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file17"
mkfifo "#file17"
cat "#file14" | head -15 > "#file16" &
cat ${IN} > "#file17" &
cat ${IN} > "#file18" &
cat "#file17" | xargs file > "#file19" &
cat "#file18" | xargs file > "#file20" &
cat "#file19" | grep "shell script" > "#file21" &
cat "#file20" | grep "shell script" > "#file22" &
cat "#file21" | cut -d: -f1 > "#file23" &
cat "#file22" | cut -d: -f1 > "#file24" &
cat "#file23" | xargs wc -l > "#file25" &
cat "#file24" | xargs wc -l > "#file26" &
cat "#file25" | sort -n > "#file27" &
cat "#file26" | sort -n > "#file28" &
sort -m --parallel=2 -n "#file27" "#file28" > "#file14" &
cat "#file16" > /tmp/distr_output/0 &
wait
rm -f "#file21"
rm -f "#file28"
rm -f "#file20"
rm -f "#file19"
rm -f "#file26"
rm -f "#file27"
rm -f "#file18"
rm -f "#file16"
rm -f "#file22"
rm -f "#file25"
rm -f "#file14"
rm -f "#file23"
rm -f "#file24"
rm -f "#file17"
rm -rf "/dev/shm/dish"