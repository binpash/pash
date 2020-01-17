rm -rf /tmp/distr_output
mkdir -p /tmp/distr_output
mkdir -p /dev/shm/dish
rm -f "#file14"
mkfifo "#file14"
rm -f "#file12"
mkfifo "#file12"
rm -f "#file17"
mkfifo "#file17"
rm -f "#file18"
mkfifo "#file18"
rm -f "#file19"
mkfifo "#file19"
rm -f "#file13"
mkfifo "#file13"
rm -f "#file11"
mkfifo "#file11"
rm -f "#file15"
mkfifo "#file15"
rm -f "#file16"
mkfifo "#file16"
cat "${IN}" > "#file11" &
cat "${IN}" > "#file12" &
cat "#file11" | tr -cs A-Za-z "\n" > "#file13" &
cat "#file12" | tr -cs A-Za-z "\n" > "#file14" &
cat "#file13" | tr A-Z a-z > "#file15" &
cat "#file14" | tr A-Z a-z > "#file16" &
alt_bigram_aux_reduce "#file17" "#file18" "#file19" &
cat "#file15" | alt_bigrams_aux  > "#file17" &
cat "#file16" | alt_bigrams_aux  > "#file18" &
cat "#file19" > /tmp/distr_output/0 &
wait
rm -f "#file14"
rm -f "#file12"
rm -f "#file17"
rm -f "#file18"
rm -f "#file19"
rm -f "#file13"
rm -f "#file11"
rm -f "#file15"
rm -f "#file16"
rm -rf "/dev/shm/dish"