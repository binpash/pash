mkfifo $t{0,1...}
curl $base/$y > $t0 &
cat $t0 | split $t1 $t2 &
cat $t1 | xargs -n1 curl -s > $t3 &
cat $t2 | xargs -n1 curl -s > $t4 &
...
cat $t9 | sort -rn > $t11 &
cat $t10 | sort -rn > $t12 &
cat $t11 | eager > $t13 &
cat $t12 | eager > $t14 &
sort -m -rn $t13 $t14 > $t15 &
cat $t15 | head -n1 > $out1 &
wait $!





mkfifo $t{0,1...}
curl $base/$y > $t0 &
cat $t0 | split $t1 $t2 &
cat $t1 | xargs -n1 curl -s > $t3 &
cat $t2 | xargs -n1 curl -s > $t4 &
cat $t3 | gunzip  > $t5 &
cat $t4 | gunzip  > $t6 &
cat $t5 | cut -c 89-92 > $t7 &
cat $t6 | cut -c 89-92 > $t8 &
cat $t7 | grep -v 999 > $t9 &
cat $t8 | grep -v 999 > $t10 &
cat $t9 | sort -rn > $t11 &
cat $t10 | sort -rn > $t12 &
cat $t11 | eager > $t13 &
cat $t12 | eager > $t14 &
sort -m -rn $t13 $t14 > $t15 &
cat $t15 | head -n1 > $out1 &
wait $!
ps --ppid $$ | awk '{print $1}' | grep -E '[0-9]' | xargs -n 1 kill -SIGPIPE
rm $t{0,1...}
