mkdir -p /tmp/dish_output/
rm -f "#file16"
mkfifo "#file16"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file19"
mkfifo "#file19"
cat $IN > "#file11" &
cat $IN > "#file12" &
cat $IN > "#file13" &
cat $IN > "#file14" &
cat "#file11" | tr A-Z a-z > "#file15" &
cat "#file12" | tr A-Z a-z > "#file16" &
cat "#file13" | tr A-Z a-z > "#file17" &
cat "#file14" | tr A-Z a-z > "#file18" &
cat "#file15" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file19" &
cat "#file16" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file20" &
cat "#file17" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file21" &
cat "#file18" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file22" &
cat "#file19" > /tmp/dish_output//0 &
cat "#file20" > /tmp/dish_output//1 &
cat "#file21" > /tmp/dish_output//2 &
cat "#file22" > /tmp/dish_output//3 &
for job in `jobs -p` 
do 
 echo $job
 wait $job 
done
rm -f "#file16"
rm -f "#file18"
rm -f "#file21"
rm -f "#file20"
rm -f "#file15"
rm -f "#file12"
rm -f "#file13"
rm -f "#file17"
rm -f "#file22"
rm -f "#file14"
rm -f "#file11"
rm -f "#file19"