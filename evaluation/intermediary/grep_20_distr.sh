rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file43"
mkfifo "#file43"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file45"
mkfifo "#file45"
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
cat ${IN} > "#file35" &
cat ${IN} > "#file36" &
cat ${IN} > "#file37" &
cat ${IN} > "#file38" &
cat ${IN} > "#file39" &
cat ${IN} > "#file40" &
cat ${IN} > "#file41" &
cat ${IN} > "#file42" &
cat ${IN} > "#file43" &
cat ${IN} > "#file44" &
cat "#file25" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file45" &
cat "#file26" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file46" &
cat "#file27" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file47" &
cat "#file28" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file48" &
cat "#file29" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file49" &
cat "#file30" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file50" &
cat "#file31" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file51" &
cat "#file32" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file52" &
cat "#file33" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file53" &
cat "#file34" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file54" &
cat "#file35" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file55" &
cat "#file36" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file56" &
cat "#file37" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file57" &
cat "#file38" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file58" &
cat "#file39" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file59" &
cat "#file40" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file60" &
cat "#file41" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file61" &
cat "#file42" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file62" &
cat "#file43" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file63" &
cat "#file44" | grep "[a-zA-Z0-9]\+@[a-zA-Z0-9]\+\.[a-z]\{2,\}" > "#file64" &
cat "#file45" > 1/0 &
cat "#file46" > 1/1 &
cat "#file47" > 1/2 &
cat "#file48" > 1/3 &
cat "#file49" > 1/4 &
cat "#file50" > 1/5 &
cat "#file51" > 1/6 &
cat "#file52" > 1/7 &
cat "#file53" > 1/8 &
cat "#file54" > 1/9 &
cat "#file55" > 1/10 &
cat "#file56" > 1/11 &
cat "#file57" > 1/12 &
cat "#file58" > 1/13 &
cat "#file59" > 1/14 &
cat "#file60" > 1/15 &
cat "#file61" > 1/16 &
cat "#file62" > 1/17 &
cat "#file63" > 1/18 &
cat "#file64" > 1/19 &
wait
rm -f "#file43"
rm -f "#file26"
rm -f "#file35"
rm -f "#file40"
rm -f "#file53"
rm -f "#file56"
rm -f "#file33"
rm -f "#file25"
rm -f "#file39"
rm -f "#file41"
rm -f "#file54"
rm -f "#file37"
rm -f "#file36"
rm -f "#file55"
rm -f "#file29"
rm -f "#file38"
rm -f "#file34"
rm -f "#file28"
rm -f "#file46"
rm -f "#file44"
rm -f "#file50"
rm -f "#file32"
rm -f "#file57"
rm -f "#file48"
rm -f "#file51"
rm -f "#file61"
rm -f "#file62"
rm -f "#file58"
rm -f "#file42"
rm -f "#file60"
rm -f "#file30"
rm -f "#file27"
rm -f "#file59"
rm -f "#file47"
rm -f "#file63"
rm -f "#file49"
rm -f "#file31"
rm -f "#file52"
rm -f "#file64"
rm -f "#file45"
rm -rf "/dev/shm/dish"