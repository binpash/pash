rm -rf /tmp/dish_output/
mkdir -p /tmp/dish_output/
mkdir -p /dev/shm/dish
rm -f "#file5"
mkfifo "#file5"
rm -f "#file7"
mkfifo "#file7"
rm -f "#file2"
mkfifo "#file2"
cat $IN > "#file2" &
cat "#file2" | tr A-Z a-z > "#file5" &
cat "#file5" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file7" &
cat "#file7" > /tmp/dish_output//0 &
wait
rm -f "#file5"
rm -f "#file7"
rm -f "#file2"
rm -rf "/dev/shm/dish"