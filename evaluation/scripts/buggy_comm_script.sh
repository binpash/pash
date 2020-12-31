mkfifo s1 s2 s3 s4

IN=../evaluation/scripts/input/1M.txt

cat $IN | sort | tee s4 >s3 &
comm -23 <(cat $IN | sort) s3 > s1 &
comm -23 <(cat $IN | sort) s4 > s2 &
cat s1 s2

rm s1 s2 s3 s4

