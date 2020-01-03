rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file73"
mkfifo "#file73"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file82"
mkfifo "#file82"
rm -f "#file84"
mkfifo "#file84"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file81"
mkfifo "#file81"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file79"
mkfifo "#file79"
rm -f "#file85"
mkfifo "#file85"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file77"
mkfifo "#file77"
rm -f "#file70"
mkfifo "#file70"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file65"
mkfifo "#file65"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file69"
mkfifo "#file69"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file78"
mkfifo "#file78"
rm -f "#file86"
mkfifo "#file86"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file67"
mkfifo "#file67"
rm -f "#file80"
mkfifo "#file80"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file83"
mkfifo "#file83"
cat "#file22" | uniq  > "#file24" &
cat "#file24" | comm -23 - ${dict} > "#file26" &
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
cat "#file27" | col -bx > "#file37" &
cat "#file28" | col -bx > "#file38" &
cat "#file29" | col -bx > "#file39" &
cat "#file30" | col -bx > "#file40" &
cat "#file31" | col -bx > "#file41" &
cat "#file32" | col -bx > "#file42" &
cat "#file33" | col -bx > "#file43" &
cat "#file34" | col -bx > "#file44" &
cat "#file35" | col -bx > "#file45" &
cat "#file36" | col -bx > "#file46" &
cat "#file37" | tr -cs A-Za-z "\n" > "#file47" &
cat "#file38" | tr -cs A-Za-z "\n" > "#file48" &
cat "#file39" | tr -cs A-Za-z "\n" > "#file49" &
cat "#file40" | tr -cs A-Za-z "\n" > "#file50" &
cat "#file41" | tr -cs A-Za-z "\n" > "#file51" &
cat "#file42" | tr -cs A-Za-z "\n" > "#file52" &
cat "#file43" | tr -cs A-Za-z "\n" > "#file53" &
cat "#file44" | tr -cs A-Za-z "\n" > "#file54" &
cat "#file45" | tr -cs A-Za-z "\n" > "#file55" &
cat "#file46" | tr -cs A-Za-z "\n" > "#file56" &
cat "#file47" | tr A-Z a-z > "#file57" &
cat "#file48" | tr A-Z a-z > "#file58" &
cat "#file49" | tr A-Z a-z > "#file59" &
cat "#file50" | tr A-Z a-z > "#file60" &
cat "#file51" | tr A-Z a-z > "#file61" &
cat "#file52" | tr A-Z a-z > "#file62" &
cat "#file53" | tr A-Z a-z > "#file63" &
cat "#file54" | tr A-Z a-z > "#file64" &
cat "#file55" | tr A-Z a-z > "#file65" &
cat "#file56" | tr A-Z a-z > "#file66" &
cat "#file57" | tr -d "[:punct:]" > "#file67" &
cat "#file58" | tr -d "[:punct:]" > "#file68" &
cat "#file59" | tr -d "[:punct:]" > "#file69" &
cat "#file60" | tr -d "[:punct:]" > "#file70" &
cat "#file61" | tr -d "[:punct:]" > "#file71" &
cat "#file62" | tr -d "[:punct:]" > "#file72" &
cat "#file63" | tr -d "[:punct:]" > "#file73" &
cat "#file64" | tr -d "[:punct:]" > "#file74" &
cat "#file65" | tr -d "[:punct:]" > "#file75" &
cat "#file66" | tr -d "[:punct:]" > "#file76" &
cat "#file67" | sort  > "#file77" &
cat "#file68" | sort  > "#file78" &
cat "#file69" | sort  > "#file79" &
cat "#file70" | sort  > "#file80" &
cat "#file71" | sort  > "#file81" &
cat "#file72" | sort  > "#file82" &
cat "#file73" | sort  > "#file83" &
cat "#file74" | sort  > "#file84" &
cat "#file75" | sort  > "#file85" &
cat "#file76" | sort  > "#file86" &
sort -m --parallel=10 "#file77" "#file78" "#file79" "#file80" "#file81" "#file82" "#file83" "#file84" "#file85" "#file86" > "#file22" &
cat "#file26" > 1/0 &
wait
rm -f "#file73"
rm -f "#file29"
rm -f "#file49"
rm -f "#file50"
rm -f "#file82"
rm -f "#file84"
rm -f "#file40"
rm -f "#file66"
rm -f "#file58"
rm -f "#file74"
rm -f "#file75"
rm -f "#file46"
rm -f "#file39"
rm -f "#file51"
rm -f "#file36"
rm -f "#file48"
rm -f "#file52"
rm -f "#file33"
rm -f "#file57"
rm -f "#file68"
rm -f "#file55"
rm -f "#file32"
rm -f "#file53"
rm -f "#file30"
rm -f "#file71"
rm -f "#file24"
rm -f "#file81"
rm -f "#file42"
rm -f "#file72"
rm -f "#file79"
rm -f "#file85"
rm -f "#file26"
rm -f "#file47"
rm -f "#file77"
rm -f "#file70"
rm -f "#file22"
rm -f "#file31"
rm -f "#file34"
rm -f "#file65"
rm -f "#file38"
rm -f "#file76"
rm -f "#file69"
rm -f "#file35"
rm -f "#file64"
rm -f "#file56"
rm -f "#file59"
rm -f "#file41"
rm -f "#file45"
rm -f "#file61"
rm -f "#file27"
rm -f "#file43"
rm -f "#file62"
rm -f "#file28"
rm -f "#file54"
rm -f "#file60"
rm -f "#file78"
rm -f "#file86"
rm -f "#file37"
rm -f "#file44"
rm -f "#file67"
rm -f "#file80"
rm -f "#file63"
rm -f "#file83"
rm -rf "/dev/shm/dish"