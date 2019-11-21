mkdir -p /tmp/dish_output/
rm -f "#file37"
mkfifo "#file37"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file78"
mkfifo "#file78"
rm -f "#file70"
mkfifo "#file70"
rm -f "#file83"
mkfifo "#file83"
rm -f "#file84"
mkfifo "#file84"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file46"
mkfifo "#file46"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file2"
mkfifo "#file2"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file65"
mkfifo "#file65"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file80"
mkfifo "#file80"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file85"
mkfifo "#file85"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file67"
mkfifo "#file67"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file86"
mkfifo "#file86"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file87"
mkfifo "#file87"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file52"
mkfifo "#file52"
cat manpage.1 > "#file2" &
cat "#file57" "#file58" "#file59" "#file60" "#file61" | uniq  > "#file15" &
cat "#file2" | split_file 100 "#file18" > "#file24" &
cat "#file24" | split_file 100 "#file19" > "#file26" &
cat "#file26" | split_file 100 "#file20" > "#file28" &
cat "#file28" | split_file 100 "#file21" > "#file22" &
cat "#file18" | groff -t -e -mandoc -Tascii > "#file31" &
cat "#file19" | groff -t -e -mandoc -Tascii > "#file32" &
cat "#file20" | groff -t -e -mandoc -Tascii > "#file33" &
cat "#file21" | groff -t -e -mandoc -Tascii > "#file34" &
cat "#file22" | groff -t -e -mandoc -Tascii > "#file35" &
cat "#file31" | col -bx > "#file36" &
cat "#file32" | col -bx > "#file37" &
cat "#file33" | col -bx > "#file38" &
cat "#file34" | col -bx > "#file39" &
cat "#file35" | col -bx > "#file40" &
cat "#file36" | tr A-Z a-z > "#file41" &
cat "#file37" | tr A-Z a-z > "#file42" &
cat "#file38" | tr A-Z a-z > "#file43" &
cat "#file39" | tr A-Z a-z > "#file44" &
cat "#file40" | tr A-Z a-z > "#file45" &
cat "#file41" | tr -d "[:punct:]" > "#file46" &
cat "#file42" | tr -d "[:punct:]" > "#file47" &
cat "#file43" | tr -d "[:punct:]" > "#file48" &
cat "#file44" | tr -d "[:punct:]" > "#file49" &
cat "#file45" | tr -d "[:punct:]" > "#file50" &
cat "#file46" | sort  > "#file51" &
cat "#file47" | sort  > "#file52" &
cat "#file48" | sort  > "#file53" &
cat "#file49" | sort  > "#file54" &
cat "#file50" | sort  > "#file55" &
sort -m "#file51" "#file52" "#file53" "#file54" "#file55" > "#file13" &
cat "#file13" | split_file 100 "#file57" > "#file63" &
cat "#file63" | split_file 100 "#file58" > "#file65" &
cat "#file65" | split_file 100 "#file59" > "#file67" &
cat "#file67" | split_file 100 "#file60" > "#file61" &
cat "#file15" | split_file 100 "#file70" > "#file76" &
cat "#file76" | split_file 100 "#file71" > "#file78" &
cat "#file78" | split_file 100 "#file72" > "#file80" &
cat "#file80" | split_file 100 "#file73" > "#file74" &
cat "#file70" | comm -13 $dict - > "#file83" &
cat "#file71" | comm -13 $dict - > "#file84" &
cat "#file72" | comm -13 $dict - > "#file85" &
cat "#file73" | comm -13 $dict - > "#file86" &
cat "#file74" | comm -13 $dict - > "#file87" &
cat "#file83" > /tmp/dish_output//0
cat "#file84" > /tmp/dish_output//1
cat "#file85" > /tmp/dish_output//2
cat "#file86" > /tmp/dish_output//3
cat "#file87" > /tmp/dish_output//4
rm -f "#file37"
rm -f "#file45"
rm -f "#file19"
rm -f "#file47"
rm -f "#file78"
rm -f "#file70"
rm -f "#file83"
rm -f "#file84"
rm -f "#file21"
rm -f "#file71"
rm -f "#file73"
rm -f "#file35"
rm -f "#file46"
rm -f "#file15"
rm -f "#file13"
rm -f "#file24"
rm -f "#file58"
rm -f "#file2"
rm -f "#file31"
rm -f "#file36"
rm -f "#file43"
rm -f "#file44"
rm -f "#file34"
rm -f "#file49"
rm -f "#file50"
rm -f "#file65"
rm -f "#file61"
rm -f "#file76"
rm -f "#file38"
rm -f "#file80"
rm -f "#file39"
rm -f "#file18"
rm -f "#file72"
rm -f "#file85"
rm -f "#file26"
rm -f "#file59"
rm -f "#file67"
rm -f "#file60"
rm -f "#file53"
rm -f "#file63"
rm -f "#file74"
rm -f "#file51"
rm -f "#file86"
rm -f "#file54"
rm -f "#file20"
rm -f "#file40"
rm -f "#file42"
rm -f "#file57"
rm -f "#file22"
rm -f "#file55"
rm -f "#file87"
rm -f "#file41"
rm -f "#file33"
rm -f "#file48"
rm -f "#file32"
rm -f "#file28"
rm -f "#file52"