rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file10"
mkfifo "#file10"
rm -f "#file8"
mkfifo "#file8"
rm -f "#file9"
mkfifo "#file9"
rm -f "#file7"
mkfifo "#file7"
cat ${IN} > "#file7" &
cat ${IN} > "#file8" &
cat "#file7" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file9" &
cat "#file8" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file10" &
cat "#file9" > 1/0 &
cat "#file10" > 1/1 &
wait
rm -f "#file10"
rm -f "#file8"
rm -f "#file9"
rm -f "#file7"
rm -rf "/dev/shm/dish"