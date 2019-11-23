rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file30"
mkfifo "#file30"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file24"
mkfifo "#file24"
cat "#file14" | uniq  > "#file16" &
cat "#file16" | comm -23 - ${dict} > "#file18" &
cat ${IN} > "#file19" &
cat ${IN} > "#file20" &
cat "#file19" | col -bx > "#file21" &
cat "#file20" | col -bx > "#file22" &
cat "#file21" | tr -cs A-Za-z "\n" > "#file23" &
cat "#file22" | tr -cs A-Za-z "\n" > "#file24" &
cat "#file23" | tr A-Z a-z > "#file25" &
cat "#file24" | tr A-Z a-z > "#file26" &
cat "#file25" | tr -d "[:punct:]" > "#file27" &
cat "#file26" | tr -d "[:punct:]" > "#file28" &
cat "#file27" | sort  > "#file29" &
cat "#file28" | sort  > "#file30" &
sort -m --parallel=2 "#file29" "#file30" > "#file14" &
cat "#file18" > 1/0 &
wait
rm -f "#file30"
rm -f "#file22"
rm -f "#file26"
rm -f "#file14"
rm -f "#file29"
rm -f "#file27"
rm -f "#file19"
rm -f "#file28"
rm -f "#file20"
rm -f "#file25"
rm -f "#file21"
rm -f "#file18"
rm -f "#file23"
rm -f "#file16"
rm -f "#file24"
rm -rf "/dev/shm/dish"