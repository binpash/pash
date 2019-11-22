rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file29"
mkfifo "#file29"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file31"
mkfifo "#file31"
cat ${IN} > "#file15" &
cat ${IN} > "#file16" &
cat ${IN} > "#file17" &
cat ${IN} > "#file18" &
cat ${IN} > "#file19" &
cat ${IN} > "#file20" &
cat ${IN} > "#file21" &
cat ${IN} > "#file22" &
cat ${IN} > "#file23" &
cat ${IN} > "#file24" &
cat "#file15" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file25" &
cat "#file16" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file26" &
cat "#file17" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file27" &
cat "#file18" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file28" &
cat "#file19" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file29" &
cat "#file20" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file30" &
cat "#file21" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file31" &
cat "#file22" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file32" &
cat "#file23" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file33" &
cat "#file24" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file34" &
cat "#file25" > /tmp/distr_output/0 &
cat "#file26" > /tmp/distr_output/1 &
cat "#file27" > /tmp/distr_output/2 &
cat "#file28" > /tmp/distr_output/3 &
cat "#file29" > /tmp/distr_output/4 &
cat "#file30" > /tmp/distr_output/5 &
cat "#file31" > /tmp/distr_output/6 &
cat "#file32" > /tmp/distr_output/7 &
cat "#file33" > /tmp/distr_output/8 &
cat "#file34" > /tmp/distr_output/9 &
wait
rm -f "#file29"
rm -f "#file28"
rm -f "#file20"
rm -f "#file16"
rm -f "#file19"
rm -f "#file17"
rm -f "#file21"
rm -f "#file23"
rm -f "#file25"
rm -f "#file30"
rm -f "#file26"
rm -f "#file32"
rm -f "#file33"
rm -f "#file27"
rm -f "#file22"
rm -f "#file24"
rm -f "#file18"
rm -f "#file15"
rm -f "#file34"
rm -f "#file31"
rm -rf "/dev/shm/dish"