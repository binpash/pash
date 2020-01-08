rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file14"
mkfifo "#file14"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file10"
mkfifo "#file10"
cat "#file10" | head -n1 > "#file12" &
cat "${IN_DIR}/p4.out_2_00" > "#file13" &
cat "${IN_DIR}/p4.out_2_01" > "#file14" &
cat "#file13" | cut -c 89-92 > "#file15" &
cat "#file14" | cut -c 89-92 > "#file16" &
cat "#file15" | grep -v 999 > "#file17" &
cat "#file16" | grep -v 999 > "#file18" &
sort -m --parallel=2 -rn "#file19" "#file20" > "#file10" &
cat "#file17" | sort -rn > "#file19" &
cat "#file18" | sort -rn > "#file20" &
cat "#file12" > /tmp/distr_output/0 &
wait
rm -f "#file14"
rm -f "#file15"
rm -f "#file18"
rm -f "#file20"
rm -f "#file16"
rm -f "#file19"
rm -f "#file12"
rm -f "#file13"
rm -f "#file17"
rm -f "#file10"
rm -rf "/dev/shm/dish"