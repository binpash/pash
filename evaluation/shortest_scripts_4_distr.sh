rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file33"
mkfifo "#file33"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file32"
mkfifo "#file32"
cat "#file16" | head -15 > "#file18" &
cat ${IN} > "#file19" &
cat ${IN} > "#file20" &
cat ${IN} > "#file21" &
cat ${IN} > "#file22" &
cat "#file19" | xargs file > "#file23" &
cat "#file20" | xargs file > "#file24" &
cat "#file21" | xargs file > "#file25" &
cat "#file22" | xargs file > "#file26" &
cat "#file23" | grep "shell script" > "#file27" &
cat "#file24" | grep "shell script" > "#file28" &
cat "#file25" | grep "shell script" > "#file29" &
cat "#file26" | grep "shell script" > "#file30" &
cat "#file27" | cut -d: -f1 > "#file31" &
cat "#file28" | cut -d: -f1 > "#file32" &
cat "#file29" | cut -d: -f1 > "#file33" &
cat "#file30" | cut -d: -f1 > "#file34" &
cat "#file31" | xargs wc -l > "#file35" &
cat "#file32" | xargs wc -l > "#file36" &
cat "#file33" | xargs wc -l > "#file37" &
cat "#file34" | xargs wc -l > "#file38" &
cat "#file35" | sort -n > "#file39" &
cat "#file36" | sort -n > "#file40" &
cat "#file37" | sort -n > "#file41" &
cat "#file38" | sort -n > "#file42" &
sort -m --parallel=4 -n "#file39" "#file40" "#file41" "#file42" > "#file16" &
cat "#file18" > /tmp/distr_output/0 &
wait
rm -f "#file33"
rm -f "#file19"
rm -f "#file38"
rm -f "#file30"
rm -f "#file23"
rm -f "#file39"
rm -f "#file27"
rm -f "#file31"
rm -f "#file20"
rm -f "#file28"
rm -f "#file40"
rm -f "#file26"
rm -f "#file18"
rm -f "#file21"
rm -f "#file37"
rm -f "#file24"
rm -f "#file16"
rm -f "#file35"
rm -f "#file29"
rm -f "#file42"
rm -f "#file36"
rm -f "#file22"
rm -f "#file41"
rm -f "#file25"
rm -f "#file34"
rm -f "#file32"
rm -rf "/dev/shm/dish"