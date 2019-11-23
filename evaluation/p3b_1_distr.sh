rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file5"
mkfifo "#file5"
rm -f "#file7"
mkfifo "#file7"
rm -f "#file2"
mkfifo "#file2"
cat "${IN_DIR}/p3a.txt" > "#file2" &
cat "#file2" | xargs -n1 curl -s > "#file5" &
cat "#file5" | gunzip  > "#file7" &
cat "#file7" > /tmp/distr_output/0 &
wait
rm -f "#file5"
rm -f "#file7"
rm -f "#file2"
rm -rf "/dev/shm/dish"