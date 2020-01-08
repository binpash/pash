rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file9"
mkfifo "#file9"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file2"
mkfifo "#file2"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file17"
mkfifo "#file17"
cat "${IN_DIR}/p4.out" > "#file2" &
cat "#file9" | head -n1 > "#file11" &
cat "#file2" | tee >( head -n 10000 > "/dev/shm/dish/#file12"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file12" > "#file12") | (tail -n +10001 > "#file13"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file12" | cut -c 89-92 > "#file16" &
cat "#file13" | cut -c 89-92 > "#file17" &
cat "#file16" | grep -v 999 > "#file18" &
cat "#file17" | grep -v 999 > "#file19" &
sort -m --parallel=2 -rn "#file20" "#file21" > "#file9" &
cat "#file18" | sort -rn > "#file20" &
cat "#file19" | sort -rn > "#file21" &
cat "#file11" > /tmp/distr_output/0 &
wait
rm -f "#file9"
rm -f "#file19"
rm -f "#file2"
rm -f "#file18"
rm -f "#file16"
rm -f "#file13"
rm -f "#file12"
rm -f "#file11"
rm -f "#file21"
rm -f "#file20"
rm -f "#file17"
rm -rf "/dev/shm/dish"