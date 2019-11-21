mkdir -p /tmp/dish_output/
rm -f "#file18"
mkfifo "#file18"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file76"
mkfifo "#file76"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file7"
mkfifo "#file7"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file51"
mkfifo "#file51"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file63"
mkfifo "#file63"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file77"
mkfifo "#file77"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file70"
mkfifo "#file70"
rm -f "#file57"
mkfifo "#file57"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file2"
mkfifo "#file2"
rm -f "#file5"
mkfifo "#file5"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file50"
mkfifo "#file50"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file48"
mkfifo "#file48"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file42"
mkfifo "#file42"
cat $IN > "#file2" &
cat "#file2" | groff -t -e -mandoc -Tascii > "#file5" &
cat "#file5" | col -bx > "#file7" &
cat "#file47" "#file48" "#file49" "#file50" "#file51" | uniq  > "#file15" &
{ head -n 100 "#file7" > "#file18" ; cat "#file7" > "#file24" ; } &
{ head -n 100 "#file24" > "#file19" ; cat "#file24" > "#file26" ; } &
{ head -n 100 "#file26" > "#file20" ; cat "#file26" > "#file28" ; } &
{ head -n 100 "#file28" > "#file21" ; cat "#file28" > "#file22" ; } &
cat "#file18" | tr A-Z a-z > "#file31" &
cat "#file19" | tr A-Z a-z > "#file32" &
cat "#file20" | tr A-Z a-z > "#file33" &
cat "#file21" | tr A-Z a-z > "#file34" &
cat "#file22" | tr A-Z a-z > "#file35" &
cat "#file31" | tr -d "[:punct:]" > "#file36" &
cat "#file32" | tr -d "[:punct:]" > "#file37" &
cat "#file33" | tr -d "[:punct:]" > "#file38" &
cat "#file34" | tr -d "[:punct:]" > "#file39" &
cat "#file35" | tr -d "[:punct:]" > "#file40" &
cat "#file36" | sort  > "#file41" &
cat "#file37" | sort  > "#file42" &
cat "#file38" | sort  > "#file43" &
cat "#file39" | sort  > "#file44" &
cat "#file40" | sort  > "#file45" &
sort -m "#file41" "#file42" "#file43" "#file44" "#file45" > "#file13" &
{ head -n 100 "#file13" > "#file47" ; cat "#file13" > "#file53" ; } &
{ head -n 100 "#file53" > "#file48" ; cat "#file53" > "#file55" ; } &
{ head -n 100 "#file55" > "#file49" ; cat "#file55" > "#file57" ; } &
{ head -n 100 "#file57" > "#file50" ; cat "#file57" > "#file51" ; } &
{ head -n 100 "#file15" > "#file60" ; cat "#file15" > "#file66" ; } &
{ head -n 100 "#file66" > "#file61" ; cat "#file66" > "#file68" ; } &
{ head -n 100 "#file68" > "#file62" ; cat "#file68" > "#file70" ; } &
{ head -n 100 "#file70" > "#file63" ; cat "#file70" > "#file64" ; } &
cat "#file60" | comm -13 $dict - > "#file73" &
cat "#file61" | comm -13 $dict - > "#file74" &
cat "#file62" | comm -13 $dict - > "#file75" &
cat "#file63" | comm -13 $dict - > "#file76" &
cat "#file64" | comm -13 $dict - > "#file77" &
cat "#file73" > /tmp/dish_output//0 &
cat "#file74" > /tmp/dish_output//1 &
cat "#file75" > /tmp/dish_output//2 &
cat "#file76" > /tmp/dish_output//3 &
cat "#file77" > /tmp/dish_output//4 &
for job in `jobs -p` 
do 
 echo $job
 wait $job 
done
rm -f "#file18"
rm -f "#file37"
rm -f "#file75"
rm -f "#file76"
rm -f "#file49"
rm -f "#file7"
rm -f "#file43"
rm -f "#file68"
rm -f "#file35"
rm -f "#file34"
rm -f "#file20"
rm -f "#file36"
rm -f "#file66"
rm -f "#file51"
rm -f "#file28"
rm -f "#file63"
rm -f "#file31"
rm -f "#file64"
rm -f "#file55"
rm -f "#file77"
rm -f "#file33"
rm -f "#file21"
rm -f "#file47"
rm -f "#file44"
rm -f "#file60"
rm -f "#file32"
rm -f "#file39"
rm -f "#file45"
rm -f "#file70"
rm -f "#file57"
rm -f "#file26"
rm -f "#file2"
rm -f "#file5"
rm -f "#file19"
rm -f "#file53"
rm -f "#file22"
rm -f "#file50"
rm -f "#file15"
rm -f "#file73"
rm -f "#file62"
rm -f "#file13"
rm -f "#file24"
rm -f "#file61"
rm -f "#file38"
rm -f "#file48"
rm -f "#file41"
rm -f "#file74"
rm -f "#file40"
rm -f "#file42"