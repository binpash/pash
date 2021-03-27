#!/bin/bash

## No newline before EOF bug
echo -n "popo" > /tmp/in
IN=/tmp/in

cat $IN $IN | grep "popopopo" > /tmp/seq.out

rm -f s1 s2
mkfifo s1 s2

cat $IN | grep "popopopo" > s1 &
cat $IN | grep "popopopo" > s2 &
cat s1 s2 > /tmp/buggy.out

rm -f s1 s2

diff /tmp/seq.out /tmp/buggy.out