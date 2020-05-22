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
sort -m -rn $t11 $t12 > $t13 &
cat $t13 | head -n1 > $out1 &
