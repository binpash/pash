rm -rf 1
mkdir -p 1
mkdir -p /dev/shm/dish
rm -f "#file101"
mkfifo "#file101"
rm -f "#file119"
mkfifo "#file119"
rm -f "#file138"
mkfifo "#file138"
rm -f "#file113"
mkfifo "#file113"
rm -f "#file97"
mkfifo "#file97"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file91"
mkfifo "#file91"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file108"
mkfifo "#file108"
rm -f "#file136"
mkfifo "#file136"
rm -f "#file148"
mkfifo "#file148"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file154"
mkfifo "#file154"
rm -f "#file100"
mkfifo "#file100"
rm -f "#file104"
mkfifo "#file104"
rm -f "#file152"
mkfifo "#file152"
rm -f "#file143"
mkfifo "#file143"
rm -f "#file153"
mkfifo "#file153"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file110"
mkfifo "#file110"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file115"
mkfifo "#file115"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file81"
mkfifo "#file81"
rm -f "#file128"
mkfifo "#file128"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file120"
mkfifo "#file120"
rm -f "#file145"
mkfifo "#file145"
rm -f "#file127"
mkfifo "#file127"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file132"
mkfifo "#file132"
rm -f "#file139"
mkfifo "#file139"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file107"
mkfifo "#file107"
rm -f "#file149"
mkfifo "#file149"
rm -f "#file147"
mkfifo "#file147"
rm -f "#file106"
mkfifo "#file106"
rm -f "#file144"
mkfifo "#file144"
rm -f "#file125"
mkfifo "#file125"
rm -f "#file129"
mkfifo "#file129"
rm -f "#file114"
mkfifo "#file114"
rm -f "#file122"
mkfifo "#file122"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file67"
mkfifo "#file67"
rm -f "#file94"
mkfifo "#file94"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file78"
mkfifo "#file78"
rm -f "#file90"
mkfifo "#file90"
rm -f "#file151"
mkfifo "#file151"
rm -f "#file105"
mkfifo "#file105"
rm -f "#file133"
mkfifo "#file133"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file80"
mkfifo "#file80"
rm -f "#file69"
mkfifo "#file69"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file123"
mkfifo "#file123"
rm -f "#file156"
mkfifo "#file156"
rm -f "#file142"
mkfifo "#file142"
rm -f "#file92"
mkfifo "#file92"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file117"
mkfifo "#file117"
rm -f "#file140"
mkfifo "#file140"
rm -f "#file102"
mkfifo "#file102"
rm -f "#file96"
mkfifo "#file96"
rm -f "#file84"
mkfifo "#file84"
rm -f "#file118"
mkfifo "#file118"
rm -f "#file150"
mkfifo "#file150"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file141"
mkfifo "#file141"
rm -f "#file130"
mkfifo "#file130"
rm -f "#file87"
mkfifo "#file87"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file116"
mkfifo "#file116"
rm -f "#file134"
mkfifo "#file134"
rm -f "#file121"
mkfifo "#file121"
rm -f "#file70"
mkfifo "#file70"
rm -f "#file155"
mkfifo "#file155"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file89"
mkfifo "#file89"
rm -f "#file103"
mkfifo "#file103"
rm -f "#file88"
mkfifo "#file88"
rm -f "#file126"
mkfifo "#file126"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file93"
mkfifo "#file93"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file146"
mkfifo "#file146"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file135"
mkfifo "#file135"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file77"
mkfifo "#file77"
rm -f "#file124"
mkfifo "#file124"
rm -f "#file131"
mkfifo "#file131"
rm -f "#file86"
mkfifo "#file86"
rm -f "#file109"
mkfifo "#file109"
rm -f "#file111"
mkfifo "#file111"
rm -f "#file79"
mkfifo "#file79"
rm -f "#file82"
mkfifo "#file82"
rm -f "#file83"
mkfifo "#file83"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file98"
mkfifo "#file98"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file99"
mkfifo "#file99"
rm -f "#file112"
mkfifo "#file112"
rm -f "#file65"
mkfifo "#file65"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file85"
mkfifo "#file85"
rm -f "#file95"
mkfifo "#file95"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file137"
mkfifo "#file137"
cat "#file32" | uniq  > "#file34" &
cat "#file34" | comm -23 - ${dict} > "#file36" &
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
cat ${IN} > "#file55" &
cat ${IN} > "#file56" &
cat "#file37" | col -bx > "#file57" &
cat "#file38" | col -bx > "#file58" &
cat "#file39" | col -bx > "#file59" &
cat "#file40" | col -bx > "#file60" &
cat "#file41" | col -bx > "#file61" &
cat "#file42" | col -bx > "#file62" &
cat "#file43" | col -bx > "#file63" &
cat "#file44" | col -bx > "#file64" &
cat "#file45" | col -bx > "#file65" &
cat "#file46" | col -bx > "#file66" &
cat "#file47" | col -bx > "#file67" &
cat "#file48" | col -bx > "#file68" &
cat "#file49" | col -bx > "#file69" &
cat "#file50" | col -bx > "#file70" &
cat "#file51" | col -bx > "#file71" &
cat "#file52" | col -bx > "#file72" &
cat "#file53" | col -bx > "#file73" &
cat "#file54" | col -bx > "#file74" &
cat "#file55" | col -bx > "#file75" &
cat "#file56" | col -bx > "#file76" &
cat "#file57" | tr -cs A-Za-z "\n" > "#file77" &
cat "#file58" | tr -cs A-Za-z "\n" > "#file78" &
cat "#file59" | tr -cs A-Za-z "\n" > "#file79" &
cat "#file60" | tr -cs A-Za-z "\n" > "#file80" &
cat "#file61" | tr -cs A-Za-z "\n" > "#file81" &
cat "#file62" | tr -cs A-Za-z "\n" > "#file82" &
cat "#file63" | tr -cs A-Za-z "\n" > "#file83" &
cat "#file64" | tr -cs A-Za-z "\n" > "#file84" &
cat "#file65" | tr -cs A-Za-z "\n" > "#file85" &
cat "#file66" | tr -cs A-Za-z "\n" > "#file86" &
cat "#file67" | tr -cs A-Za-z "\n" > "#file87" &
cat "#file68" | tr -cs A-Za-z "\n" > "#file88" &
cat "#file69" | tr -cs A-Za-z "\n" > "#file89" &
cat "#file70" | tr -cs A-Za-z "\n" > "#file90" &
cat "#file71" | tr -cs A-Za-z "\n" > "#file91" &
cat "#file72" | tr -cs A-Za-z "\n" > "#file92" &
cat "#file73" | tr -cs A-Za-z "\n" > "#file93" &
cat "#file74" | tr -cs A-Za-z "\n" > "#file94" &
cat "#file75" | tr -cs A-Za-z "\n" > "#file95" &
cat "#file76" | tr -cs A-Za-z "\n" > "#file96" &
cat "#file77" | tr A-Z a-z > "#file97" &
cat "#file78" | tr A-Z a-z > "#file98" &
cat "#file79" | tr A-Z a-z > "#file99" &
cat "#file80" | tr A-Z a-z > "#file100" &
cat "#file81" | tr A-Z a-z > "#file101" &
cat "#file82" | tr A-Z a-z > "#file102" &
cat "#file83" | tr A-Z a-z > "#file103" &
cat "#file84" | tr A-Z a-z > "#file104" &
cat "#file85" | tr A-Z a-z > "#file105" &
cat "#file86" | tr A-Z a-z > "#file106" &
cat "#file87" | tr A-Z a-z > "#file107" &
cat "#file88" | tr A-Z a-z > "#file108" &
cat "#file89" | tr A-Z a-z > "#file109" &
cat "#file90" | tr A-Z a-z > "#file110" &
cat "#file91" | tr A-Z a-z > "#file111" &
cat "#file92" | tr A-Z a-z > "#file112" &
cat "#file93" | tr A-Z a-z > "#file113" &
cat "#file94" | tr A-Z a-z > "#file114" &
cat "#file95" | tr A-Z a-z > "#file115" &
cat "#file96" | tr A-Z a-z > "#file116" &
cat "#file97" | tr -d "[:punct:]" > "#file117" &
cat "#file98" | tr -d "[:punct:]" > "#file118" &
cat "#file99" | tr -d "[:punct:]" > "#file119" &
cat "#file100" | tr -d "[:punct:]" > "#file120" &
cat "#file101" | tr -d "[:punct:]" > "#file121" &
cat "#file102" | tr -d "[:punct:]" > "#file122" &
cat "#file103" | tr -d "[:punct:]" > "#file123" &
cat "#file104" | tr -d "[:punct:]" > "#file124" &
cat "#file105" | tr -d "[:punct:]" > "#file125" &
cat "#file106" | tr -d "[:punct:]" > "#file126" &
cat "#file107" | tr -d "[:punct:]" > "#file127" &
cat "#file108" | tr -d "[:punct:]" > "#file128" &
cat "#file109" | tr -d "[:punct:]" > "#file129" &
cat "#file110" | tr -d "[:punct:]" > "#file130" &
cat "#file111" | tr -d "[:punct:]" > "#file131" &
cat "#file112" | tr -d "[:punct:]" > "#file132" &
cat "#file113" | tr -d "[:punct:]" > "#file133" &
cat "#file114" | tr -d "[:punct:]" > "#file134" &
cat "#file115" | tr -d "[:punct:]" > "#file135" &
cat "#file116" | tr -d "[:punct:]" > "#file136" &
cat "#file117" | sort  > "#file137" &
cat "#file118" | sort  > "#file138" &
cat "#file119" | sort  > "#file139" &
cat "#file120" | sort  > "#file140" &
cat "#file121" | sort  > "#file141" &
cat "#file122" | sort  > "#file142" &
cat "#file123" | sort  > "#file143" &
cat "#file124" | sort  > "#file144" &
cat "#file125" | sort  > "#file145" &
cat "#file126" | sort  > "#file146" &
cat "#file127" | sort  > "#file147" &
cat "#file128" | sort  > "#file148" &
cat "#file129" | sort  > "#file149" &
cat "#file130" | sort  > "#file150" &
cat "#file131" | sort  > "#file151" &
cat "#file132" | sort  > "#file152" &
cat "#file133" | sort  > "#file153" &
cat "#file134" | sort  > "#file154" &
cat "#file135" | sort  > "#file155" &
cat "#file136" | sort  > "#file156" &
sort -m --parallel=20 "#file137" "#file138" "#file139" "#file140" "#file141" "#file142" "#file143" "#file144" "#file145" "#file146" "#file147" "#file148" "#file149" "#file150" "#file151" "#file152" "#file153" "#file154" "#file155" "#file156" > "#file32" &
cat "#file36" > 1/0 &
wait
rm -f "#file101"
rm -f "#file119"
rm -f "#file138"
rm -f "#file113"
rm -f "#file97"
rm -f "#file54"
rm -f "#file91"
rm -f "#file56"
rm -f "#file75"
rm -f "#file76"
rm -f "#file57"
rm -f "#file108"
rm -f "#file136"
rm -f "#file148"
rm -f "#file74"
rm -f "#file154"
rm -f "#file100"
rm -f "#file104"
rm -f "#file152"
rm -f "#file143"
rm -f "#file153"
rm -f "#file41"
rm -f "#file110"
rm -f "#file60"
rm -f "#file50"
rm -f "#file115"
rm -f "#file61"
rm -f "#file81"
rm -f "#file128"
rm -f "#file55"
rm -f "#file36"
rm -f "#file64"
rm -f "#file48"
rm -f "#file34"
rm -f "#file120"
rm -f "#file145"
rm -f "#file127"
rm -f "#file45"
rm -f "#file132"
rm -f "#file139"
rm -f "#file38"
rm -f "#file107"
rm -f "#file149"
rm -f "#file147"
rm -f "#file106"
rm -f "#file144"
rm -f "#file125"
rm -f "#file129"
rm -f "#file114"
rm -f "#file122"
rm -f "#file32"
rm -f "#file43"
rm -f "#file67"
rm -f "#file94"
rm -f "#file52"
rm -f "#file78"
rm -f "#file90"
rm -f "#file151"
rm -f "#file105"
rm -f "#file133"
rm -f "#file51"
rm -f "#file80"
rm -f "#file69"
rm -f "#file47"
rm -f "#file123"
rm -f "#file156"
rm -f "#file142"
rm -f "#file92"
rm -f "#file59"
rm -f "#file62"
rm -f "#file40"
rm -f "#file117"
rm -f "#file140"
rm -f "#file102"
rm -f "#file96"
rm -f "#file84"
rm -f "#file118"
rm -f "#file150"
rm -f "#file68"
rm -f "#file141"
rm -f "#file130"
rm -f "#file87"
rm -f "#file72"
rm -f "#file116"
rm -f "#file134"
rm -f "#file121"
rm -f "#file70"
rm -f "#file155"
rm -f "#file37"
rm -f "#file89"
rm -f "#file103"
rm -f "#file88"
rm -f "#file126"
rm -f "#file71"
rm -f "#file93"
rm -f "#file46"
rm -f "#file146"
rm -f "#file63"
rm -f "#file42"
rm -f "#file135"
rm -f "#file49"
rm -f "#file77"
rm -f "#file124"
rm -f "#file131"
rm -f "#file86"
rm -f "#file109"
rm -f "#file111"
rm -f "#file79"
rm -f "#file82"
rm -f "#file83"
rm -f "#file44"
rm -f "#file98"
rm -f "#file58"
rm -f "#file53"
rm -f "#file99"
rm -f "#file112"
rm -f "#file65"
rm -f "#file73"
rm -f "#file66"
rm -f "#file85"
rm -f "#file95"
rm -f "#file39"
rm -f "#file137"
rm -rf "/dev/shm/dish"