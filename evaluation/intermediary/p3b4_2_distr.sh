rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file19"
mkfifo "#file19"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file17"
mkfifo "#file17"
cat "#file14" | head -n1 > "#file16" &
cat "${IN_DIR}/p3a_2_1.txt" > "#file17" &
cat "${IN_DIR}/p3a_2_2.txt" > "#file18" &
cat "#file17" | xargs -n1 curl -s > "#file19" &
cat "#file18" | xargs -n1 curl -s > "#file20" &
cat "#file19" | gunzip  > "#file21" &
cat "#file20" | gunzip  > "#file22" &
cat "#file21" | cut -c 89-92 > "#file23" &
cat "#file22" | cut -c 89-92 > "#file24" &
cat "#file23" | grep -v 999 > "#file25" &
cat "#file24" | grep -v 999 > "#file26" &
cat "#file25" | sort -rn > "#file27" &
cat "#file26" | sort -rn > "#file28" &
sort -m --parallel=2 -rn "#file27" "#file28" > "#file14" &
cat "#file16" > /tmp/distr_output/0 &
wait
rm -f "#file19"
rm -f "#file22"
rm -f "#file18"
rm -f "#file14"
rm -f "#file24"
rm -f "#file25"
rm -f "#file26"
rm -f "#file27"
rm -f "#file28"
rm -f "#file21"
rm -f "#file16"
rm -f "#file23"
rm -f "#file20"
rm -f "#file17"
rm -rf "/dev/shm/dish"