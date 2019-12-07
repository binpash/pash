rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file26"
mkfifo "#file26"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file32"
mkfifo "#file32"
cat ${IN} > "#file17" &
cat ${IN} > "#file18" &
cat ${IN} > "#file19" &
cat ${IN} > "#file20" &
cat ${IN} > "#file21" &
cat ${IN} > "#file22" &
cat ${IN} > "#file23" &
cat ${IN} > "#file24" &
cat ${IN} > "#file25" &
cat ${IN} > "#file26" &
cat "#file17" | tr -cs A-Za-z "\n" > "#file27" &
cat "#file18" | tr -cs A-Za-z "\n" > "#file28" &
cat "#file19" | tr -cs A-Za-z "\n" > "#file29" &
cat "#file20" | tr -cs A-Za-z "\n" > "#file30" &
cat "#file21" | tr -cs A-Za-z "\n" > "#file31" &
cat "#file22" | tr -cs A-Za-z "\n" > "#file32" &
cat "#file23" | tr -cs A-Za-z "\n" > "#file33" &
cat "#file24" | tr -cs A-Za-z "\n" > "#file34" &
cat "#file25" | tr -cs A-Za-z "\n" > "#file35" &
cat "#file26" | tr -cs A-Za-z "\n" > "#file36" &
cat "#file27" | tr A-Z a-z > "#file37" &
cat "#file28" | tr A-Z a-z > "#file38" &
cat "#file29" | tr A-Z a-z > "#file39" &
cat "#file30" | tr A-Z a-z > "#file40" &
cat "#file31" | tr A-Z a-z > "#file41" &
cat "#file32" | tr A-Z a-z > "#file42" &
cat "#file33" | tr A-Z a-z > "#file43" &
cat "#file34" | tr A-Z a-z > "#file44" &
cat "#file35" | tr A-Z a-z > "#file45" &
cat "#file36" | tr A-Z a-z > "#file46" &
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
rm -f "#file26"
rm -f "#file34"
rm -f "#file42"
rm -f "#file45"
rm -f "#file28"
rm -f "#file35"
rm -f "#file22"
rm -f "#file46"
rm -f "#file21"
rm -f "#file43"
rm -f "#file27"
rm -f "#file41"
rm -f "#file37"
rm -f "#file19"
rm -f "#file40"
rm -f "#file33"
rm -f "#file30"
rm -f "#file44"
rm -f "#file24"
rm -f "#file23"
rm -f "#file25"
rm -f "#file31"
rm -f "#file29"
rm -f "#file36"
rm -f "#file17"
rm -f "#file18"
rm -f "#file20"
rm -f "#file39"
rm -f "#file38"
rm -f "#file32"
rm -rf "/dev/shm/dish"