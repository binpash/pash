rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file37"
mkfifo "#file37"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file20"
mkfifo "#file20"
cat "#file16" | uniq  > "#file18" &
cat "#file18" | comm -23 - ${dict} > "#file20" &
cat ${IN} > "#file21" &
cat ${IN} > "#file22" &
cat ${IN} > "#file23" &
cat ${IN} > "#file24" &
cat "#file21" | col -bx > "#file25" &
cat "#file22" | col -bx > "#file26" &
cat "#file23" | col -bx > "#file27" &
cat "#file24" | col -bx > "#file28" &
cat "#file25" | tr -cs A-Za-z "\n" > "#file29" &
cat "#file26" | tr -cs A-Za-z "\n" > "#file30" &
cat "#file27" | tr -cs A-Za-z "\n" > "#file31" &
cat "#file28" | tr -cs A-Za-z "\n" > "#file32" &
cat "#file29" | tr A-Z a-z > "#file33" &
cat "#file30" | tr A-Z a-z > "#file34" &
cat "#file31" | tr A-Z a-z > "#file35" &
cat "#file32" | tr A-Z a-z > "#file36" &
cat "#file33" | tr -d "[:punct:]" > "#file37" &
cat "#file34" | tr -d "[:punct:]" > "#file38" &
cat "#file35" | tr -d "[:punct:]" > "#file39" &
cat "#file36" | tr -d "[:punct:]" > "#file40" &
cat "#file37" | sort  > "#file41" &
cat "#file38" | sort  > "#file42" &
cat "#file39" | sort  > "#file43" &
cat "#file40" | sort  > "#file44" &
sort -m --parallel=4 "#file41" "#file42" "#file43" "#file44" > "#file16" &
cat "#file20" > 1/0 &
wait
rm -f "#file37"
rm -f "#file38"
rm -f "#file34"
rm -f "#file29"
rm -f "#file21"
rm -f "#file31"
rm -f "#file22"
rm -f "#file25"
rm -f "#file35"
rm -f "#file44"
rm -f "#file16"
rm -f "#file36"
rm -f "#file28"
rm -f "#file43"
rm -f "#file33"
rm -f "#file41"
rm -f "#file27"
rm -f "#file26"
rm -f "#file32"
rm -f "#file18"
rm -f "#file30"
rm -f "#file40"
rm -f "#file42"
rm -f "#file39"
rm -f "#file23"
rm -f "#file24"
rm -f "#file20"
rm -rf "/dev/shm/dish"