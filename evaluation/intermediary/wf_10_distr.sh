rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file20"
mkfifo "#file20"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file60"
mkfifo "#file60"
cat "#file18" | uniq -c > "#file20" &
cat "#file20" | sort -rn > "#file22" &
cat "#file22" | sed ${1}q > "#file24" &
cat ${IN} > "#file25" &
cat ${IN} > "#file26" &
cat ${IN} > "#file27" &
cat ${IN} > "#file28" &
cat ${IN} > "#file29" &
cat ${IN} > "#file30" &
cat ${IN} > "#file31" &
cat ${IN} > "#file32" &
cat ${IN} > "#file33" &
cat ${IN} > "#file34" &
cat "#file25" | tr -cs A-Za-z"\n" > "#file35" &
cat "#file26" | tr -cs A-Za-z"\n" > "#file36" &
cat "#file27" | tr -cs A-Za-z"\n" > "#file37" &
cat "#file28" | tr -cs A-Za-z"\n" > "#file38" &
cat "#file29" | tr -cs A-Za-z"\n" > "#file39" &
cat "#file30" | tr -cs A-Za-z"\n" > "#file40" &
cat "#file31" | tr -cs A-Za-z"\n" > "#file41" &
cat "#file32" | tr -cs A-Za-z"\n" > "#file42" &
cat "#file33" | tr -cs A-Za-z"\n" > "#file43" &
cat "#file34" | tr -cs A-Za-z"\n" > "#file44" &
cat "#file35" | tr A-Z a-z > "#file45" &
cat "#file36" | tr A-Z a-z > "#file46" &
cat "#file37" | tr A-Z a-z > "#file47" &
cat "#file38" | tr A-Z a-z > "#file48" &
cat "#file39" | tr A-Z a-z > "#file49" &
cat "#file40" | tr A-Z a-z > "#file50" &
cat "#file41" | tr A-Z a-z > "#file51" &
cat "#file42" | tr A-Z a-z > "#file52" &
cat "#file43" | tr A-Z a-z > "#file53" &
cat "#file44" | tr A-Z a-z > "#file54" &
cat "#file45" | sort  > "#file55" &
cat "#file46" | sort  > "#file56" &
cat "#file47" | sort  > "#file57" &
cat "#file48" | sort  > "#file58" &
cat "#file49" | sort  > "#file59" &
cat "#file50" | sort  > "#file60" &
cat "#file51" | sort  > "#file61" &
cat "#file52" | sort  > "#file62" &
cat "#file53" | sort  > "#file63" &
cat "#file54" | sort  > "#file64" &
sort -m --parallel=10 "#file55" "#file56" "#file57" "#file58" "#file59" "#file60" "#file61" "#file62" "#file63" "#file64" > "#file18" &
cat "#file24" > 1/0 &
wait
rm -f "#file20"
rm -f "#file59"
rm -f "#file35"
rm -f "#file47"
rm -f "#file57"
rm -f "#file38"
rm -f "#file34"
rm -f "#file52"
rm -f "#file49"
rm -f "#file22"
rm -f "#file33"
rm -f "#file62"
rm -f "#file43"
rm -f "#file63"
rm -f "#file25"
rm -f "#file37"
rm -f "#file61"
rm -f "#file30"
rm -f "#file51"
rm -f "#file42"
rm -f "#file45"
rm -f "#file18"
rm -f "#file55"
rm -f "#file32"
rm -f "#file58"
rm -f "#file28"
rm -f "#file46"
rm -f "#file64"
rm -f "#file56"
rm -f "#file31"
rm -f "#file39"
rm -f "#file54"
rm -f "#file29"
rm -f "#file48"
rm -f "#file36"
rm -f "#file27"
rm -f "#file24"
rm -f "#file26"
rm -f "#file41"
rm -f "#file40"
rm -f "#file44"
rm -f "#file53"
rm -f "#file50"
rm -f "#file60"
rm -rf "/dev/shm/dish"