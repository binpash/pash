#!/bin/bash
echo “Enter Any Number”
read n
i=1
c=1
while [ $i -le $n ]; do
	i=$(expr $i + 1)
	r=$(expr $n % $i)
	if [ $r -eq 0 ]; then
		c=$(expr $c + 1)
	fi
done
if [ $c -eq 2 ]; then
	echo “Prime”
else
	echo “Not Prime”
fi
