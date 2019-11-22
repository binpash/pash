rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file2"
mkfifo "#file2"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file5"
mkfifo "#file5"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file7"
mkfifo "#file7"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file9"
mkfifo "#file9"
cat ${IN} > "#file2" &
cat "#file2" | col -bx > "#file5" &
cat "#file5" | tr -cs A-Za-z "\n" > "#file7" &
cat "#file7" | tr A-Z a-z > "#file9" &
cat "#file9" | tr -d "[:punct:]" > "#file11" &
cat "#file11" | sort  > "#file13" &
cat "#file13" | uniq  > "#file15" &
cat "#file15" | comm -23 - ${dict} > "#file17" &
cat "#file17" > /tmp/distr_output/0 &
wait
rm -f "#file2"
rm -f "#file17"
rm -f "#file11"
rm -f "#file5"
rm -f "#file13"
rm -f "#file7"
rm -f "#file15"
rm -f "#file9"
rm -rf "/dev/shm/dish"