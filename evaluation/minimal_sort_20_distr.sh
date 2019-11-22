rm -rf /tmp/dish_output/
mkdir -p /tmp/dish_output/
mkdir -p /dev/shm/dish
rm -f "#file70"
mkfifo "#file70"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file67"
mkfifo "#file67"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file82"
mkfifo "#file82"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file69"
mkfifo "#file69"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file77"
mkfifo "#file77"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file80"
mkfifo "#file80"
rm -f "#file83"
mkfifo "#file83"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file78"
mkfifo "#file78"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file65"
mkfifo "#file65"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file79"
mkfifo "#file79"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file86"
mkfifo "#file86"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file84"
mkfifo "#file84"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file85"
mkfifo "#file85"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file81"
mkfifo "#file81"
rm -f "#file57"
mkfifo "#file57"
cat $IN > "#file27" &
cat $IN > "#file28" &
cat $IN > "#file29" &
cat $IN > "#file30" &
cat $IN > "#file31" &
cat $IN > "#file32" &
cat $IN > "#file33" &
cat $IN > "#file34" &
cat $IN > "#file35" &
cat $IN > "#file36" &
cat $IN > "#file37" &
cat $IN > "#file38" &
cat $IN > "#file39" &
cat $IN > "#file40" &
cat $IN > "#file41" &
cat $IN > "#file42" &
cat $IN > "#file43" &
cat $IN > "#file44" &
cat $IN > "#file45" &
cat $IN > "#file46" &
cat "#file27" | tr A-Z a-z > "#file47" &
cat "#file28" | tr A-Z a-z > "#file48" &
cat "#file29" | tr A-Z a-z > "#file49" &
cat "#file30" | tr A-Z a-z > "#file50" &
cat "#file31" | tr A-Z a-z > "#file51" &
cat "#file32" | tr A-Z a-z > "#file52" &
cat "#file33" | tr A-Z a-z > "#file53" &
cat "#file34" | tr A-Z a-z > "#file54" &
cat "#file35" | tr A-Z a-z > "#file55" &
cat "#file36" | tr A-Z a-z > "#file56" &
cat "#file37" | tr A-Z a-z > "#file57" &
cat "#file38" | tr A-Z a-z > "#file58" &
cat "#file39" | tr A-Z a-z > "#file59" &
cat "#file40" | tr A-Z a-z > "#file60" &
cat "#file41" | tr A-Z a-z > "#file61" &
cat "#file42" | tr A-Z a-z > "#file62" &
cat "#file43" | tr A-Z a-z > "#file63" &
cat "#file44" | tr A-Z a-z > "#file64" &
cat "#file45" | tr A-Z a-z > "#file65" &
cat "#file46" | tr A-Z a-z > "#file66" &
cat "#file47" | sort  > "#file67" &
cat "#file48" | sort  > "#file68" &
cat "#file49" | sort  > "#file69" &
cat "#file50" | sort  > "#file70" &
cat "#file51" | sort  > "#file71" &
cat "#file52" | sort  > "#file72" &
cat "#file53" | sort  > "#file73" &
cat "#file54" | sort  > "#file74" &
cat "#file55" | sort  > "#file75" &
cat "#file56" | sort  > "#file76" &
cat "#file57" | sort  > "#file77" &
cat "#file58" | sort  > "#file78" &
cat "#file59" | sort  > "#file79" &
cat "#file60" | sort  > "#file80" &
cat "#file61" | sort  > "#file81" &
cat "#file62" | sort  > "#file82" &
cat "#file63" | sort  > "#file83" &
cat "#file64" | sort  > "#file84" &
cat "#file65" | sort  > "#file85" &
cat "#file66" | sort  > "#file86" &
sort -m --parallel=20 "#file67" "#file68" "#file69" "#file70" "#file71" "#file72" "#file73" "#file74" "#file75" "#file76" "#file77" "#file78" "#file79" "#file80" "#file81" "#file82" "#file83" "#file84" "#file85" "#file86" > "#file26" &
cat "#file26" > /tmp/dish_output//0 &
wait
rm -f "#file70"
rm -f "#file59"
rm -f "#file60"
rm -f "#file62"
rm -f "#file75"
rm -f "#file61"
rm -f "#file67"
rm -f "#file30"
rm -f "#file74"
rm -f "#file29"
rm -f "#file27"
rm -f "#file39"
rm -f "#file33"
rm -f "#file82"
rm -f "#file66"
rm -f "#file35"
rm -f "#file69"
rm -f "#file34"
rm -f "#file51"
rm -f "#file56"
rm -f "#file64"
rm -f "#file77"
rm -f "#file28"
rm -f "#file37"
rm -f "#file40"
rm -f "#file49"
rm -f "#file80"
rm -f "#file83"
rm -f "#file26"
rm -f "#file55"
rm -f "#file50"
rm -f "#file78"
rm -f "#file48"
rm -f "#file63"
rm -f "#file42"
rm -f "#file47"
rm -f "#file32"
rm -f "#file65"
rm -f "#file58"
rm -f "#file79"
rm -f "#file41"
rm -f "#file44"
rm -f "#file45"
rm -f "#file86"
rm -f "#file31"
rm -f "#file68"
rm -f "#file36"
rm -f "#file71"
rm -f "#file43"
rm -f "#file53"
rm -f "#file73"
rm -f "#file84"
rm -f "#file76"
rm -f "#file54"
rm -f "#file72"
rm -f "#file52"
rm -f "#file85"
rm -f "#file46"
rm -f "#file38"
rm -f "#file81"
rm -f "#file57"
rm -rf "/dev/shm/dish"