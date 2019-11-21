mkdir -p /tmp/dish_output/
rm -f "#file46"
mkfifo "#file46"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file43"
mkfifo "#file43"
rm -f "#file33"
mkfifo "#file33"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file37"
mkfifo "#file37"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file35"
mkfifo "#file35"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file38"
mkfifo "#file38"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file44"
mkfifo "#file44"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file21"
mkfifo "#file21"
cat $IN > "#file17" &
cat $IN > "#file18" &
cat $IN > "#file19" &
cat $IN > "#file20" &
cat $IN > "#file21" &
cat $IN > "#file22" &
cat $IN > "#file23" &
cat $IN > "#file24" &
cat $IN > "#file25" &
cat $IN > "#file26" &
cat "#file17" | tr A-Z a-z > "#file27" &
cat "#file18" | tr A-Z a-z > "#file28" &
cat "#file19" | tr A-Z a-z > "#file29" &
cat "#file20" | tr A-Z a-z > "#file30" &
cat "#file21" | tr A-Z a-z > "#file31" &
cat "#file22" | tr A-Z a-z > "#file32" &
cat "#file23" | tr A-Z a-z > "#file33" &
cat "#file24" | tr A-Z a-z > "#file34" &
cat "#file25" | tr A-Z a-z > "#file35" &
cat "#file26" | tr A-Z a-z > "#file36" &
cat "#file27" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file37" &
cat "#file28" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file38" &
cat "#file29" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file39" &
cat "#file30" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file40" &
cat "#file31" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file41" &
cat "#file32" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file42" &
cat "#file33" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file43" &
cat "#file34" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file44" &
cat "#file35" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file45" &
cat "#file36" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file46" &
cat "#file37" > /tmp/dish_output//0 &
cat "#file38" > /tmp/dish_output//1 &
cat "#file39" > /tmp/dish_output//2 &
cat "#file40" > /tmp/dish_output//3 &
cat "#file41" > /tmp/dish_output//4 &
cat "#file42" > /tmp/dish_output//5 &
cat "#file43" > /tmp/dish_output//6 &
cat "#file44" > /tmp/dish_output//7 &
cat "#file45" > /tmp/dish_output//8 &
cat "#file46" > /tmp/dish_output//9 &
for job in `jobs -p` 
do 
 echo $job
 wait $job 
done
rm -f "#file46"
rm -f "#file32"
rm -f "#file41"
rm -f "#file43"
rm -f "#file33"
rm -f "#file30"
rm -f "#file34"
rm -f "#file37"
rm -f "#file23"
rm -f "#file29"
rm -f "#file45"
rm -f "#file20"
rm -f "#file25"
rm -f "#file18"
rm -f "#file35"
rm -f "#file26"
rm -f "#file22"
rm -f "#file38"
rm -f "#file27"
rm -f "#file28"
rm -f "#file17"
rm -f "#file31"
rm -f "#file36"
rm -f "#file39"
rm -f "#file19"
rm -f "#file42"
rm -f "#file40"
rm -f "#file44"
rm -f "#file24"
rm -f "#file21"