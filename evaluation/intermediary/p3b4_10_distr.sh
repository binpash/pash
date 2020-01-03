rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file30"
mkfifo "#file30"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file79"
mkfifo "#file79"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file70"
mkfifo "#file70"
rm -f "#file81"
mkfifo "#file81"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file80"
mkfifo "#file80"
rm -f "#file84"
mkfifo "#file84"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file69"
mkfifo "#file69"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file65"
mkfifo "#file65"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file78"
mkfifo "#file78"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file67"
mkfifo "#file67"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file77"
mkfifo "#file77"
rm -f "#file83"
mkfifo "#file83"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file82"
mkfifo "#file82"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file42"
mkfifo "#file42"
cat "#file22" | head -n1 > "#file24" &
cat "${IN_DIR}/p3a_10_1.txt" > "#file25" &
cat "${IN_DIR}/p3a_10_2.txt" > "#file26" &
cat "${IN_DIR}/p3a_10_3.txt" > "#file27" &
cat "${IN_DIR}/p3a_10_4.txt" > "#file28" &
cat "${IN_DIR}/p3a_10_5.txt" > "#file29" &
cat "${IN_DIR}/p3a_10_6.txt" > "#file30" &
cat "${IN_DIR}/p3a_10_7.txt" > "#file31" &
cat "${IN_DIR}/p3a_10_8.txt" > "#file32" &
cat "${IN_DIR}/p3a_10_9.txt" > "#file33" &
cat "${IN_DIR}/p3a_10_10.txt" > "#file34" &
cat "#file25" | xargs -n1 curl -s > "#file35" &
cat "#file26" | xargs -n1 curl -s > "#file36" &
cat "#file27" | xargs -n1 curl -s > "#file37" &
cat "#file28" | xargs -n1 curl -s > "#file38" &
cat "#file29" | xargs -n1 curl -s > "#file39" &
cat "#file30" | xargs -n1 curl -s > "#file40" &
cat "#file31" | xargs -n1 curl -s > "#file41" &
cat "#file32" | xargs -n1 curl -s > "#file42" &
cat "#file33" | xargs -n1 curl -s > "#file43" &
cat "#file34" | xargs -n1 curl -s > "#file44" &
cat "#file35" | gunzip  > "#file45" &
cat "#file36" | gunzip  > "#file46" &
cat "#file37" | gunzip  > "#file47" &
cat "#file38" | gunzip  > "#file48" &
cat "#file39" | gunzip  > "#file49" &
cat "#file40" | gunzip  > "#file50" &
cat "#file41" | gunzip  > "#file51" &
cat "#file42" | gunzip  > "#file52" &
cat "#file43" | gunzip  > "#file53" &
cat "#file44" | gunzip  > "#file54" &
cat "#file45" | cut -c 89-92 > "#file55" &
cat "#file46" | cut -c 89-92 > "#file56" &
cat "#file47" | cut -c 89-92 > "#file57" &
cat "#file48" | cut -c 89-92 > "#file58" &
cat "#file49" | cut -c 89-92 > "#file59" &
cat "#file50" | cut -c 89-92 > "#file60" &
cat "#file51" | cut -c 89-92 > "#file61" &
cat "#file52" | cut -c 89-92 > "#file62" &
cat "#file53" | cut -c 89-92 > "#file63" &
cat "#file54" | cut -c 89-92 > "#file64" &
cat "#file55" | grep -v 999 > "#file65" &
cat "#file56" | grep -v 999 > "#file66" &
cat "#file57" | grep -v 999 > "#file67" &
cat "#file58" | grep -v 999 > "#file68" &
cat "#file59" | grep -v 999 > "#file69" &
cat "#file60" | grep -v 999 > "#file70" &
cat "#file61" | grep -v 999 > "#file71" &
cat "#file62" | grep -v 999 > "#file72" &
cat "#file63" | grep -v 999 > "#file73" &
cat "#file64" | grep -v 999 > "#file74" &
cat "#file65" | sort -rn > "#file75" &
cat "#file66" | sort -rn > "#file76" &
cat "#file67" | sort -rn > "#file77" &
cat "#file68" | sort -rn > "#file78" &
cat "#file69" | sort -rn > "#file79" &
cat "#file70" | sort -rn > "#file80" &
cat "#file71" | sort -rn > "#file81" &
cat "#file72" | sort -rn > "#file82" &
cat "#file73" | sort -rn > "#file83" &
cat "#file74" | sort -rn > "#file84" &
sort -m --parallel=10 -rn "#file75" "#file76" "#file77" "#file78" "#file79" "#file80" "#file81" "#file82" "#file83" "#file84" > "#file22" &
cat "#file24" > /tmp/distr_output/0 &
wait
rm -f "#file30"
rm -f "#file57"
rm -f "#file60"
rm -f "#file34"
rm -f "#file24"
rm -f "#file79"
rm -f "#file64"
rm -f "#file54"
rm -f "#file22"
rm -f "#file71"
rm -f "#file36"
rm -f "#file31"
rm -f "#file62"
rm -f "#file26"
rm -f "#file40"
rm -f "#file73"
rm -f "#file74"
rm -f "#file70"
rm -f "#file81"
rm -f "#file43"
rm -f "#file76"
rm -f "#file45"
rm -f "#file35"
rm -f "#file38"
rm -f "#file58"
rm -f "#file59"
rm -f "#file80"
rm -f "#file84"
rm -f "#file29"
rm -f "#file41"
rm -f "#file69"
rm -f "#file27"
rm -f "#file39"
rm -f "#file65"
rm -f "#file25"
rm -f "#file56"
rm -f "#file50"
rm -f "#file66"
rm -f "#file78"
rm -f "#file72"
rm -f "#file51"
rm -f "#file67"
rm -f "#file53"
rm -f "#file55"
rm -f "#file49"
rm -f "#file75"
rm -f "#file33"
rm -f "#file77"
rm -f "#file83"
rm -f "#file68"
rm -f "#file32"
rm -f "#file61"
rm -f "#file46"
rm -f "#file37"
rm -f "#file52"
rm -f "#file82"
rm -f "#file44"
rm -f "#file48"
rm -f "#file63"
rm -f "#file47"
rm -f "#file28"
rm -f "#file42"
rm -rf "/dev/shm/dish"