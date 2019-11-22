rm -rf /tmp/dish_output/
mkdir -p /tmp/dish_output/
mkdir -p /dev/shm/dish
rm -f "#file28"
mkfifo "#file28"
rm -f "#file36"
mkfifo "#file36"
rm -f "#file47"
mkfifo "#file47"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file68"
mkfifo "#file68"
rm -f "#file72"
mkfifo "#file72"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file39"
mkfifo "#file39"
rm -f "#file34"
mkfifo "#file34"
rm -f "#file71"
mkfifo "#file71"
rm -f "#file59"
mkfifo "#file59"
rm -f "#file62"
mkfifo "#file62"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file27"
mkfifo "#file27"
rm -f "#file55"
mkfifo "#file55"
rm -f "#file75"
mkfifo "#file75"
rm -f "#file21"
mkfifo "#file21"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file10"
mkfifo "#file10"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file53"
mkfifo "#file53"
rm -f "#file40"
mkfifo "#file40"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file66"
mkfifo "#file66"
rm -f "#file49"
mkfifo "#file49"
rm -f "#file42"
mkfifo "#file42"
rm -f "#file52"
mkfifo "#file52"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file64"
mkfifo "#file64"
rm -f "#file58"
mkfifo "#file58"
rm -f "#file73"
mkfifo "#file73"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file56"
mkfifo "#file56"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file60"
mkfifo "#file60"
rm -f "#file45"
mkfifo "#file45"
rm -f "#file74"
mkfifo "#file74"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file61"
mkfifo "#file61"
rm -f "#file54"
mkfifo "#file54"
rm -f "#file41"
mkfifo "#file41"
rm -f "#file43"
mkfifo "#file43"
cat "#file26" "#file27" "#file28" "#file29" "#file30" | uniq -c > "#file12" &
cat $IN > "#file17" &
cat $IN > "#file18" &
cat "#file17" | tr -cs A-Za-z "\n" > "#file19" &
cat "#file18" | tr -cs A-Za-z "\n" > "#file20" &
cat "#file19" | tr A-Z a-z > "#file21" &
cat "#file20" | tr A-Z a-z > "#file22" &
cat "#file21" | sort  > "#file23" &
cat "#file22" | sort  > "#file24" &
sort -m --parallel=2 "#file23" "#file24" > "#file10" &
cat "#file10" | tee >( head -n 10000 > "/dev/shm/dish/#file26"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file26" > "#file26") | (tail -n +10001 > "#file32"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file32" | tee >( head -n 10000 > "/dev/shm/dish/#file27"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file27" > "#file27") | (tail -n +10001 > "#file34"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file34" | tee >( head -n 10000 > "/dev/shm/dish/#file28"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file28" > "#file28") | (tail -n +10001 > "#file36"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file36" | tee >( head -n 10000 > "/dev/shm/dish/#file29"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file29" > "#file29") | (tail -n +10001 > "#file30"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file12" | tee >( head -n 10000 > "/dev/shm/dish/#file39"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file39" > "#file39") | (tail -n +10001 > "#file45"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file45" | tee >( head -n 10000 > "/dev/shm/dish/#file40"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file40" > "#file40") | (tail -n +10001 > "#file47"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file47" | tee >( head -n 10000 > "/dev/shm/dish/#file41"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file41" > "#file41") | (tail -n +10001 > "#file49"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file49" | tee >( head -n 10000 > "/dev/shm/dish/#file42"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file42" > "#file42") | (tail -n +10001 > "#file43"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file39" | sort -rn > "#file52" &
cat "#file40" | sort -rn > "#file53" &
cat "#file41" | sort -rn > "#file54" &
cat "#file42" | sort -rn > "#file55" &
cat "#file43" | sort -rn > "#file56" &
sort -m --parallel=5 "#file52" "#file53" "#file54" "#file55" "#file56" > "#file14" &
cat "#file14" | tee >( head -n 10000 > "/dev/shm/dish/#file58"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file58" > "#file58") | (tail -n +10001 > "#file64"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file64" | tee >( head -n 10000 > "/dev/shm/dish/#file59"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file59" > "#file59") | (tail -n +10001 > "#file66"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file66" | tee >( head -n 10000 > "/dev/shm/dish/#file60"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file60" > "#file60") | (tail -n +10001 > "#file68"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file68" | tee >( head -n 10000 > "/dev/shm/dish/#file61"; dd of=/dev/null > /dev/null 2>&1 & cat "/dev/shm/dish/#file61" > "#file61") | (tail -n +10001 > "#file62"; dd of=/dev/null > /dev/null 2>&1) &
cat "#file58" | sed $Nq > "#file71" &
cat "#file59" | sed $Nq > "#file72" &
cat "#file60" | sed $Nq > "#file73" &
cat "#file61" | sed $Nq > "#file74" &
cat "#file62" | sed $Nq > "#file75" &
cat "#file71" > /tmp/dish_output//0 &
cat "#file72" > /tmp/dish_output//1 &
cat "#file73" > /tmp/dish_output//2 &
cat "#file74" > /tmp/dish_output//3 &
cat "#file75" > /tmp/dish_output//4 &
wait
rm -f "#file28"
rm -f "#file36"
rm -f "#file47"
rm -f "#file26"
rm -f "#file68"
rm -f "#file72"
rm -f "#file32"
rm -f "#file12"
rm -f "#file18"
rm -f "#file39"
rm -f "#file34"
rm -f "#file71"
rm -f "#file59"
rm -f "#file62"
rm -f "#file23"
rm -f "#file27"
rm -f "#file55"
rm -f "#file75"
rm -f "#file21"
rm -f "#file17"
rm -f "#file10"
rm -f "#file24"
rm -f "#file53"
rm -f "#file40"
rm -f "#file19"
rm -f "#file66"
rm -f "#file49"
rm -f "#file42"
rm -f "#file52"
rm -f "#file29"
rm -f "#file64"
rm -f "#file58"
rm -f "#file73"
rm -f "#file30"
rm -f "#file56"
rm -f "#file20"
rm -f "#file14"
rm -f "#file60"
rm -f "#file45"
rm -f "#file74"
rm -f "#file22"
rm -f "#file61"
rm -f "#file54"
rm -f "#file41"
rm -f "#file43"
rm -rf "/dev/shm/dish"