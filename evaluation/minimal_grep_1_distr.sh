mkdir -p /tmp/dish_output/
rm -f "#file14"
mkfifo "#file14"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file8"
mkfifo "#file8"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file9"
mkfifo "#file9"
rm -f "#file2"
mkfifo "#file2"
rm -f "#file13"
mkfifo "#file13"
cat $IN > "#file2" &
{ head -n 100 "#file2" > "#file8" ; cat "#file2" > "#file9" ; } &
cat "#file8" | tr A-Z a-z > "#file12" &
cat "#file9" | tr A-Z a-z > "#file13" &
cat "#file12" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file14" &
cat "#file13" | grep "\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4" > "#file15" &
cat "#file14" > /tmp/dish_output//0 &
cat "#file15" > /tmp/dish_output//1 &
for job in `jobs -p` 
do 
 echo $job
 wait $job 
done
rm -f "#file14"
rm -f "#file12"
rm -f "#file8"
rm -f "#file15"
rm -f "#file9"
rm -f "#file2"
rm -f "#file13"