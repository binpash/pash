rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file10"
mkfifo "#file10"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file9"
mkfifo "#file9"
cat ${IN} > "#file9" &
cat ${IN} > "#file10" &
cat ${IN} > "#file11" &
cat ${IN} > "#file12" &
cat "#file9" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file13" &
cat "#file10" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file14" &
cat "#file11" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file15" &
cat "#file12" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file16" &
cat "#file13" > 1/0 &
cat "#file14" > 1/1 &
cat "#file15" > 1/2 &
cat "#file16" > 1/3 &
wait
rm -f "#file10"
rm -f "#file12"
rm -f "#file13"
rm -f "#file14"
rm -f "#file16"
rm -f "#file15"
rm -f "#file11"
rm -f "#file9"
rm -rf "/dev/shm/dish"