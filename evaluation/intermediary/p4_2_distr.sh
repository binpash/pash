rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file13"
mkfifo "#file13"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file9"
mkfifo "#file9"
rm -f "#file10"
mkfifo "#file10"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file12"
mkfifo "#file12"
cat "${IN_DIR}/p3.out_2_00" > "#file9" &
cat "${IN_DIR}/p3.out_2_01" > "#file10" &
cat "#file9" | xargs -n1 curl -s > "#file11" &
cat "#file10" | xargs -n1 curl -s > "#file12" &
cat "#file11" | gunzip  > "#file13" &
cat "#file12" | gunzip  > "#file14" &
cat "#file13" > /tmp/distr_output/0 &
cat "#file14" > /tmp/distr_output/1 &
wait
rm -f "#file13"
rm -f "#file14"
rm -f "#file9"
rm -f "#file10"
rm -f "#file11"
rm -f "#file12"
rm -rf "/dev/shm/dish"