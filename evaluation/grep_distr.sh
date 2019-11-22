rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file4"
mkfifo "#file4"
rm -f "#file28"
mkfifo "#file28"
rm -f "#file20"
mkfifo "#file20"
rm -f "#file22"
mkfifo "#file22"
rm -f "#file14"
mkfifo "#file14"
rm -f "#file24"
mkfifo "#file24"
rm -f "#file8"
mkfifo "#file8"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file10"
mkfifo "#file10"
rm -f "#file6"
mkfifo "#file6"
rm -f "#file1"
mkfifo "#file1"
rm -f "#file30"
mkfifo "#file30"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file16"
mkfifo "#file16"
rm -f "#file2"
mkfifo "#file2"
rm -f "#file26"
mkfifo "#file26"
cat "#file1" | seq 2005 2005 > "#file2" &
cat "#file2" | sed "s;^;http://ndr.md/noaa/;" > "#file4" &
cat "#file4" | sed "s;$;/;" > "#file6" &
cat "#file6" | xargs -n 1 curl -s > "#file8" &
cat "#file8" | grep gz > "#file10" &
cat "#file10" | head -n 1 > "#file12" &
cat "#file12" | tr -s " " > "#file14" &
cat "#file14" | cut -d " " -f9 > "#file16" &
cat "#file16" | sed "s;^;http://ndr.md/noaa/2005/;" > "#file18" &
cat "#file18" | xargs -n1 curl -s > "#file20" &
cat "#file20" | gunzip  > "#file22" &
cat "#file22" | cut -c 89-92 > "#file24" &
cat "#file24" | grep -v 999 > "#file26" &
cat "#file26" | sort -rn > "#file28" &
cat "#file28" | head -n1 > "#file30" &
cat "#file30" > /tmp/distr_output/0 &
wait
rm -f "#file4"
rm -f "#file28"
rm -f "#file20"
rm -f "#file22"
rm -f "#file14"
rm -f "#file24"
rm -f "#file8"
rm -f "#file18"
rm -f "#file10"
rm -f "#file6"
rm -f "#file1"
rm -f "#file30"
rm -f "#file12"
rm -f "#file16"
rm -f "#file2"
rm -f "#file26"
rm -rf "/dev/shm/dish"