mkfifo s1 s2

cat ../evaluation/scripts/input/1M.txt > s1 &
cat ../evaluation/scripts/input/1M.txt > s2 &
(cat s1 s2; head -n 1 s2 > /dev/null) | head -n 1 &

wait

rm s1 s2
