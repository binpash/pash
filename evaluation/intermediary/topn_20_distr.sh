rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file86"
mkfifo "#file86"
rm -f "#file88"
mkfifo "#file88"
rm -f "#file102"
mkfifo "#file102"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file106"
mkfifo "#file106"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file69"
mkfifo "#file69"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file94"
mkfifo "#file94"
rm -f "#file67"
mkfifo "#file67"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file93"
mkfifo "#file93"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file84"
mkfifo "#file84"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file105"
mkfifo "#file105"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file113"
mkfifo "#file113"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file109"
mkfifo "#file109"
rm -f "#file108"
mkfifo "#file108"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file89"
mkfifo "#file89"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file87"
mkfifo "#file87"
rm -f "#file107"
mkfifo "#file107"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file78"
mkfifo "#file78"
rm -f "#file65"
mkfifo "#file65"
rm -f "#file110"
mkfifo "#file110"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file98"
mkfifo "#file98"
rm -f "#file70"
mkfifo "#file70"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file91"
mkfifo "#file91"
rm -f "#file92"
mkfifo "#file92"
rm -f "#file85"
mkfifo "#file85"
rm -f "#file90"
mkfifo "#file90"
rm -f "#file104"
mkfifo "#file104"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file99"
mkfifo "#file99"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file96"
mkfifo "#file96"
rm -f "#file97"
mkfifo "#file97"
rm -f "#file101"
mkfifo "#file101"
rm -f "#file112"
mkfifo "#file112"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file103"
mkfifo "#file103"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file114"
mkfifo "#file114"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file111"
mkfifo "#file111"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file83"
mkfifo "#file83"
rm -f "#file79"
mkfifo "#file79"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file81"
mkfifo "#file81"
rm -f "#file82"
mkfifo "#file82"
rm -f "#file77"
mkfifo "#file77"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file100"
mkfifo "#file100"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file95"
mkfifo "#file95"
rm -f "#file80"
mkfifo "#file80"
cat "#file28" | uniq -c > "#file30" &
cat "#file30" | sort -rn > "#file32" &
cat "#file32" | sed ${N}q > "#file34" &
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
cat ${IN} > "#file45" &
cat ${IN} > "#file46" &
cat ${IN} > "#file47" &
cat ${IN} > "#file48" &
cat ${IN} > "#file49" &
cat ${IN} > "#file50" &
cat ${IN} > "#file51" &
cat ${IN} > "#file52" &
cat ${IN} > "#file53" &
cat ${IN} > "#file54" &
cat "#file35" | tr -cs A-Za-z "\n" > "#file55" &
cat "#file36" | tr -cs A-Za-z "\n" > "#file56" &
cat "#file37" | tr -cs A-Za-z "\n" > "#file57" &
cat "#file38" | tr -cs A-Za-z "\n" > "#file58" &
cat "#file39" | tr -cs A-Za-z "\n" > "#file59" &
cat "#file40" | tr -cs A-Za-z "\n" > "#file60" &
cat "#file41" | tr -cs A-Za-z "\n" > "#file61" &
cat "#file42" | tr -cs A-Za-z "\n" > "#file62" &
cat "#file43" | tr -cs A-Za-z "\n" > "#file63" &
cat "#file44" | tr -cs A-Za-z "\n" > "#file64" &
cat "#file45" | tr -cs A-Za-z "\n" > "#file65" &
cat "#file46" | tr -cs A-Za-z "\n" > "#file66" &
cat "#file47" | tr -cs A-Za-z "\n" > "#file67" &
cat "#file48" | tr -cs A-Za-z "\n" > "#file68" &
cat "#file49" | tr -cs A-Za-z "\n" > "#file69" &
cat "#file50" | tr -cs A-Za-z "\n" > "#file70" &
cat "#file51" | tr -cs A-Za-z "\n" > "#file71" &
cat "#file52" | tr -cs A-Za-z "\n" > "#file72" &
cat "#file53" | tr -cs A-Za-z "\n" > "#file73" &
cat "#file54" | tr -cs A-Za-z "\n" > "#file74" &
cat "#file55" | tr A-Z a-z > "#file75" &
cat "#file56" | tr A-Z a-z > "#file76" &
cat "#file57" | tr A-Z a-z > "#file77" &
cat "#file58" | tr A-Z a-z > "#file78" &
cat "#file59" | tr A-Z a-z > "#file79" &
cat "#file60" | tr A-Z a-z > "#file80" &
cat "#file61" | tr A-Z a-z > "#file81" &
cat "#file62" | tr A-Z a-z > "#file82" &
cat "#file63" | tr A-Z a-z > "#file83" &
cat "#file64" | tr A-Z a-z > "#file84" &
cat "#file65" | tr A-Z a-z > "#file85" &
cat "#file66" | tr A-Z a-z > "#file86" &
cat "#file67" | tr A-Z a-z > "#file87" &
cat "#file68" | tr A-Z a-z > "#file88" &
cat "#file69" | tr A-Z a-z > "#file89" &
cat "#file70" | tr A-Z a-z > "#file90" &
cat "#file71" | tr A-Z a-z > "#file91" &
cat "#file72" | tr A-Z a-z > "#file92" &
cat "#file73" | tr A-Z a-z > "#file93" &
cat "#file74" | tr A-Z a-z > "#file94" &
cat "#file75" | sort  > "#file95" &
cat "#file76" | sort  > "#file96" &
cat "#file77" | sort  > "#file97" &
cat "#file78" | sort  > "#file98" &
cat "#file79" | sort  > "#file99" &
cat "#file80" | sort  > "#file100" &
cat "#file81" | sort  > "#file101" &
cat "#file82" | sort  > "#file102" &
cat "#file83" | sort  > "#file103" &
cat "#file84" | sort  > "#file104" &
cat "#file85" | sort  > "#file105" &
cat "#file86" | sort  > "#file106" &
cat "#file87" | sort  > "#file107" &
cat "#file88" | sort  > "#file108" &
cat "#file89" | sort  > "#file109" &
cat "#file90" | sort  > "#file110" &
cat "#file91" | sort  > "#file111" &
cat "#file92" | sort  > "#file112" &
cat "#file93" | sort  > "#file113" &
cat "#file94" | sort  > "#file114" &
sort -m --parallel=20 "#file95" "#file96" "#file97" "#file98" "#file99" "#file100" "#file101" "#file102" "#file103" "#file104" "#file105" "#file106" "#file107" "#file108" "#file109" "#file110" "#file111" "#file112" "#file113" "#file114" > "#file28" &
cat "#file34" > 1/0 &
wait
rm -f "#file86"
rm -f "#file88"
rm -f "#file102"
rm -f "#file45"
rm -f "#file106"
rm -f "#file63"
rm -f "#file69"
rm -f "#file62"
rm -f "#file74"
rm -f "#file41"
rm -f "#file94"
rm -f "#file67"
rm -f "#file60"
rm -f "#file61"
rm -f "#file54"
rm -f "#file53"
rm -f "#file55"
rm -f "#file76"
rm -f "#file57"
rm -f "#file93"
rm -f "#file34"
rm -f "#file84"
rm -f "#file73"
rm -f "#file43"
rm -f "#file105"
rm -f "#file71"
rm -f "#file30"
rm -f "#file113"
rm -f "#file32"
rm -f "#file109"
rm -f "#file108"
rm -f "#file64"
rm -f "#file89"
rm -f "#file50"
rm -f "#file38"
rm -f "#file58"
rm -f "#file40"
rm -f "#file87"
rm -f "#file107"
rm -f "#file68"
rm -f "#file52"
rm -f "#file78"
rm -f "#file65"
rm -f "#file110"
rm -f "#file56"
rm -f "#file98"
rm -f "#file70"
rm -f "#file72"
rm -f "#file91"
rm -f "#file92"
rm -f "#file85"
rm -f "#file90"
rm -f "#file104"
rm -f "#file49"
rm -f "#file46"
rm -f "#file99"
rm -f "#file36"
rm -f "#file96"
rm -f "#file97"
rm -f "#file101"
rm -f "#file112"
rm -f "#file35"
rm -f "#file103"
rm -f "#file51"
rm -f "#file114"
rm -f "#file37"
rm -f "#file39"
rm -f "#file111"
rm -f "#file66"
rm -f "#file28"
rm -f "#file83"
rm -f "#file79"
rm -f "#file48"
rm -f "#file81"
rm -f "#file82"
rm -f "#file77"
rm -f "#file42"
rm -f "#file44"
rm -f "#file100"
rm -f "#file47"
rm -f "#file59"
rm -f "#file75"
rm -f "#file95"
rm -f "#file80"
rm -rf "/dev/shm/dish"