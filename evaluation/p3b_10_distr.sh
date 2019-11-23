rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file21"
mkfifo "#file21"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file24"
mkfifo "#file24"
cat "${IN_DIR}/p3a_10_1.txt" > "#file17" &
cat "${IN_DIR}/p3a_10_2.txt" > "#file18" &
cat "${IN_DIR}/p3a_10_3.txt" > "#file19" &
cat "${IN_DIR}/p3a_10_4.txt" > "#file20" &
cat "${IN_DIR}/p3a_10_5.txt" > "#file21" &
cat "${IN_DIR}/p3a_10_6.txt" > "#file22" &
cat "${IN_DIR}/p3a_10_7.txt" > "#file23" &
cat "${IN_DIR}/p3a_10_8.txt" > "#file24" &
cat "${IN_DIR}/p3a_10_9.txt" > "#file25" &
cat "${IN_DIR}/p3a_10_10.txt" > "#file26" &
cat "#file17" | xargs -n1 curl -s > "#file27" &
cat "#file18" | xargs -n1 curl -s > "#file28" &
cat "#file19" | xargs -n1 curl -s > "#file29" &
cat "#file20" | xargs -n1 curl -s > "#file30" &
cat "#file21" | xargs -n1 curl -s > "#file31" &
cat "#file22" | xargs -n1 curl -s > "#file32" &
cat "#file23" | xargs -n1 curl -s > "#file33" &
cat "#file24" | xargs -n1 curl -s > "#file34" &
cat "#file25" | xargs -n1 curl -s > "#file35" &
cat "#file26" | xargs -n1 curl -s > "#file36" &
cat "#file27" | gunzip  > "#file37" &
cat "#file28" | gunzip  > "#file38" &
cat "#file29" | gunzip  > "#file39" &
cat "#file30" | gunzip  > "#file40" &
cat "#file31" | gunzip  > "#file41" &
cat "#file32" | gunzip  > "#file42" &
cat "#file33" | gunzip  > "#file43" &
cat "#file34" | gunzip  > "#file44" &
cat "#file35" | gunzip  > "#file45" &
cat "#file36" | gunzip  > "#file46" &
cat "#file37" > /tmp/distr_output/0 &
cat "#file38" > /tmp/distr_output/1 &
cat "#file39" > /tmp/distr_output/2 &
cat "#file40" > /tmp/distr_output/3 &
cat "#file41" > /tmp/distr_output/4 &
cat "#file42" > /tmp/distr_output/5 &
cat "#file43" > /tmp/distr_output/6 &
cat "#file44" > /tmp/distr_output/7 &
cat "#file45" > /tmp/distr_output/8 &
cat "#file46" > /tmp/distr_output/9 &
wait
rm -f "#file21"
rm -f "#file34"
rm -f "#file36"
rm -f "#file44"
rm -f "#file31"
rm -f "#file33"
rm -f "#file37"
rm -f "#file19"
rm -f "#file46"
rm -f "#file27"
rm -f "#file45"
rm -f "#file41"
rm -f "#file29"
rm -f "#file18"
rm -f "#file23"
rm -f "#file26"
rm -f "#file38"
rm -f "#file17"
rm -f "#file30"
rm -f "#file28"
rm -f "#file40"
rm -f "#file39"
rm -f "#file35"
rm -f "#file25"
rm -f "#file20"
rm -f "#file32"
rm -f "#file42"
rm -f "#file43"
rm -f "#file22"
rm -f "#file24"
rm -rf "/dev/shm/dish"