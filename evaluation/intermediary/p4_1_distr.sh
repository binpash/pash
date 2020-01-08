rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file9"
mkfifo "#file9"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file2"
mkfifo "#file2"
rm -f "#file8"
mkfifo "#file8"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file14"
mkfifo "#file14"
cat "${IN_DIR}/p3.out" > "#file2" &
cat "#file2" | tee >( head -n 10000 > "/dev/shm/dish/#file8"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file8" > "#file8") | (tail -n +10001 > "#file9"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file8" | xargs -n1 curl -s > "#file12" &
cat "#file9" | xargs -n1 curl -s > "#file13" &
cat "#file12" | gunzip  > "#file14" &
cat "#file13" | gunzip  > "#file15" &
cat "#file14" > /tmp/distr_output/0 &
cat "#file15" > /tmp/distr_output/1 &
wait
rm -f "#file9"
rm -f "#file15"
rm -f "#file13"
rm -f "#file2"
rm -f "#file8"
rm -f "#file12"
rm -f "#file14"
rm -rf "/dev/shm/dish"