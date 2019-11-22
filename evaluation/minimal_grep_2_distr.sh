rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file14"
mkfifo "#file14"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file9"
mkfifo "#file9"
rm -f "#file10"
mkfifo "#file10"
rm -f "#file12"
mkfifo "#file12"
cat ${IN} > "#file9" &
cat ${IN} > "#file10" &
cat "#file9" | tr A-Z a-z > "#file11" &
cat "#file10" | tr A-Z a-z > "#file12" &
cat "#file11" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file13" &
cat "#file12" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file14" &
cat "#file13" > 1/0 &
cat "#file14" > 1/1 &
wait
rm -f "#file14"
rm -f "#file13"
rm -f "#file11"
rm -f "#file9"
rm -f "#file10"
rm -f "#file12"
rm -rf "/dev/shm/dish"