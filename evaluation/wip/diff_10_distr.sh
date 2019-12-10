rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file22"
mkfifo "#file22"
rm -f "#file25"
mkfifo "#file25"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file29"
mkfifo "#file29"
rm -f "#file32"
mkfifo "#file32"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file23"
mkfifo "#file23"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file31"
mkfifo "#file31"
rm -f "#file26"
mkfifo "#file26"
rm -f "#file27"
mkfifo "#file27"
cat "#file23" "#file24" "#file25" "#file26" "#file27" "#file28" "#file29" "#file30" "#file31" "#file32" | awk 'BEGIN {srand(); OFMT="%.17f"} {print rand(), $0}' '${@}' > "#file14" &
cat "#file14" | sort -k1,1n > "#file16" &
cat "#file16" | cut -d ' ' -f2- > "#file18" &
cat "#file18" | sort  > "#file20" &
cat "#file20" | tr [:lower:] [:upper:] > "#file22" &
cat ${IN} > "#file23" &
cat ${IN} > "#file24" &
cat ${IN} > "#file25" &
cat ${IN} > "#file26" &
cat ${IN} > "#file27" &
cat ${IN} > "#file28" &
cat ${IN} > "#file29" &
cat ${IN} > "#file30" &
cat ${IN} > "#file31" &
cat ${IN} > "#file32" &
cat "#file22" > /tmp/distr_output/0 &
wait
rm -f "#file22"
rm -f "#file25"
rm -f "#file16"
rm -f "#file20"
rm -f "#file29"
rm -f "#file32"
rm -f "#file28"
rm -f "#file23"
rm -f "#file24"
rm -f "#file18"
rm -f "#file30"
rm -f "#file14"
rm -f "#file31"
rm -f "#file26"
rm -f "#file27"
rm -rf "/dev/shm/dish"