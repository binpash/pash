rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file30"
mkfifo "#file30"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file78"
mkfifo "#file78"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file67"
mkfifo "#file67"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file82"
mkfifo "#file82"
rm -f "#file69"
mkfifo "#file69"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file77"
mkfifo "#file77"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file65"
mkfifo "#file65"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file84"
mkfifo "#file84"
rm -f "#file80"
mkfifo "#file80"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file81"
mkfifo "#file81"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file70"
mkfifo "#file70"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file83"
mkfifo "#file83"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file79"
mkfifo "#file79"
cat "#file22" | head -15 > "#file24" &
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
cat "#file25" | xargs file > "#file35" &
cat "#file26" | xargs file > "#file36" &
cat "#file27" | xargs file > "#file37" &
cat "#file28" | xargs file > "#file38" &
cat "#file29" | xargs file > "#file39" &
cat "#file30" | xargs file > "#file40" &
cat "#file31" | xargs file > "#file41" &
cat "#file32" | xargs file > "#file42" &
cat "#file33" | xargs file > "#file43" &
cat "#file34" | xargs file > "#file44" &
cat "#file35" | grep "shell script" > "#file45" &
cat "#file36" | grep "shell script" > "#file46" &
cat "#file37" | grep "shell script" > "#file47" &
cat "#file38" | grep "shell script" > "#file48" &
cat "#file39" | grep "shell script" > "#file49" &
cat "#file40" | grep "shell script" > "#file50" &
cat "#file41" | grep "shell script" > "#file51" &
cat "#file42" | grep "shell script" > "#file52" &
cat "#file43" | grep "shell script" > "#file53" &
cat "#file44" | grep "shell script" > "#file54" &
cat "#file45" | cut -d: -f1 > "#file55" &
cat "#file46" | cut -d: -f1 > "#file56" &
cat "#file47" | cut -d: -f1 > "#file57" &
cat "#file48" | cut -d: -f1 > "#file58" &
cat "#file49" | cut -d: -f1 > "#file59" &
cat "#file50" | cut -d: -f1 > "#file60" &
cat "#file51" | cut -d: -f1 > "#file61" &
cat "#file52" | cut -d: -f1 > "#file62" &
cat "#file53" | cut -d: -f1 > "#file63" &
cat "#file54" | cut -d: -f1 > "#file64" &
cat "#file55" | xargs wc -l > "#file65" &
cat "#file56" | xargs wc -l > "#file66" &
cat "#file57" | xargs wc -l > "#file67" &
cat "#file58" | xargs wc -l > "#file68" &
cat "#file59" | xargs wc -l > "#file69" &
cat "#file60" | xargs wc -l > "#file70" &
cat "#file61" | xargs wc -l > "#file71" &
cat "#file62" | xargs wc -l > "#file72" &
cat "#file63" | xargs wc -l > "#file73" &
cat "#file64" | xargs wc -l > "#file74" &
cat "#file65" | sort -n > "#file75" &
cat "#file66" | sort -n > "#file76" &
cat "#file67" | sort -n > "#file77" &
cat "#file68" | sort -n > "#file78" &
cat "#file69" | sort -n > "#file79" &
cat "#file70" | sort -n > "#file80" &
cat "#file71" | sort -n > "#file81" &
cat "#file72" | sort -n > "#file82" &
cat "#file73" | sort -n > "#file83" &
cat "#file74" | sort -n > "#file84" &
sort -m --parallel=10 -n "#file75" "#file76" "#file77" "#file78" "#file79" "#file80" "#file81" "#file82" "#file83" "#file84" > "#file22" &
cat "#file24" > /tmp/distr_output/0 &
wait
rm -f "#file30"
rm -f "#file45"
rm -f "#file40"
rm -f "#file64"
rm -f "#file78"
rm -f "#file39"
rm -f "#file34"
rm -f "#file48"
rm -f "#file35"
rm -f "#file52"
rm -f "#file32"
rm -f "#file22"
rm -f "#file67"
rm -f "#file38"
rm -f "#file31"
rm -f "#file73"
rm -f "#file58"
rm -f "#file62"
rm -f "#file71"
rm -f "#file82"
rm -f "#file69"
rm -f "#file68"
rm -f "#file26"
rm -f "#file41"
rm -f "#file75"
rm -f "#file36"
rm -f "#file77"
rm -f "#file61"
rm -f "#file27"
rm -f "#file59"
rm -f "#file37"
rm -f "#file24"
rm -f "#file47"
rm -f "#file54"
rm -f "#file65"
rm -f "#file66"
rm -f "#file72"
rm -f "#file84"
rm -f "#file80"
rm -f "#file33"
rm -f "#file25"
rm -f "#file74"
rm -f "#file81"
rm -f "#file43"
rm -f "#file55"
rm -f "#file56"
rm -f "#file49"
rm -f "#file57"
rm -f "#file60"
rm -f "#file70"
rm -f "#file76"
rm -f "#file83"
rm -f "#file28"
rm -f "#file51"
rm -f "#file63"
rm -f "#file29"
rm -f "#file44"
rm -f "#file42"
rm -f "#file50"
rm -f "#file53"
rm -f "#file46"
rm -f "#file79"
rm -rf "/dev/shm/dish"