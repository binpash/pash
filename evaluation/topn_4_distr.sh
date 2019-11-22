rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file25"
mkfifo "#file25"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file22"
mkfifo "#file22"
cat "#file12" | uniq -c > "#file14" &
cat "#file14" | sort -rn > "#file16" &
cat "#file16" | sed ${N}q > "#file18" &
cat ${IN} > "#file19" &
cat ${IN} > "#file20" &
cat ${IN} > "#file21" &
cat ${IN} > "#file22" &
cat "#file19" | tr -cs A-Za-z "\n" > "#file23" &
cat "#file20" | tr -cs A-Za-z "\n" > "#file24" &
cat "#file21" | tr -cs A-Za-z "\n" > "#file25" &
cat "#file22" | tr -cs A-Za-z "\n" > "#file26" &
cat "#file23" | tr A-Z a-z > "#file27" &
cat "#file24" | tr A-Z a-z > "#file28" &
cat "#file25" | tr A-Z a-z > "#file29" &
cat "#file26" | tr A-Z a-z > "#file30" &
cat "#file27" | sort  > "#file31" &
cat "#file28" | sort  > "#file32" &
cat "#file29" | sort  > "#file33" &
cat "#file30" | sort  > "#file34" &
sort -m --parallel=4 "#file31" "#file32" "#file33" "#file34" > "#file12" &
cat "#file18" > /tmp/distr_output/0 &
wait
rm -f "#file25"
rm -f "#file19"
rm -f "#file28"
rm -f "#file32"
rm -f "#file21"
rm -f "#file14"
rm -f "#file29"
rm -f "#file12"
rm -f "#file16"
rm -f "#file18"
rm -f "#file20"
rm -f "#file23"
rm -f "#file27"
rm -f "#file30"
rm -f "#file34"
rm -f "#file31"
rm -f "#file24"
rm -f "#file33"
rm -f "#file26"
rm -f "#file22"
rm -rf "/dev/shm/dish"