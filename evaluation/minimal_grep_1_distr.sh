mkdir -p /tmp/dish_output/
rm -f "#file28"
mkfifo "#file28"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file8"
mkfifo "#file8"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file10"
mkfifo "#file10"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file2"
mkfifo "#file2"
rm -f "#file9"
mkfifo "#file9"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file27"
mkfifo "#file27"
cat $IN > "#file2" &
{ head -n 100 "#file2" > "#file8" ; cat "#file2" > "#file14" ; } &
{ head -n 100 "#file14" > "#file9" ; cat "#file14" > "#file16" ; } &
{ head -n 100 "#file16" > "#file10" ; cat "#file16" > "#file18" ; } &
{ head -n 100 "#file18" > "#file11" ; cat "#file18" > "#file12" ; } &
cat "#file8" | tr A-Z a-z > "#file21" &
cat "#file9" | tr A-Z a-z > "#file22" &
cat "#file10" | tr A-Z a-z > "#file23" &
cat "#file11" | tr A-Z a-z > "#file24" &
cat "#file12" | tr A-Z a-z > "#file25" &
cat "#file21" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file26" &
cat "#file22" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file27" &
cat "#file23" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file28" &
cat "#file24" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file29" &
cat "#file25" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file30" &
cat "#file26" > /tmp/dish_output//0 &
cat "#file27" > /tmp/dish_output//1 &
cat "#file28" > /tmp/dish_output//2 &
cat "#file29" > /tmp/dish_output//3 &
cat "#file30" > /tmp/dish_output//4 &
for job in `jobs -p` 
do 
 echo $job
 wait $job 
done
rm -f "#file28"
rm -f "#file18"
rm -f "#file22"
rm -f "#file14"
rm -f "#file8"
rm -f "#file21"
rm -f "#file11"
rm -f "#file23"
rm -f "#file12"
rm -f "#file25"
rm -f "#file10"
rm -f "#file30"
rm -f "#file24"
rm -f "#file16"
rm -f "#file29"
rm -f "#file2"
rm -f "#file9"
rm -f "#file26"
rm -f "#file27"