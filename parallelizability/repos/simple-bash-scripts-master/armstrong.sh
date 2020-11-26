#!/bin/bash
echo "Enter A Number"
read n
arm=0
temp=$n
while [ $n -ne 0 ]; do
	r=$(expr $n % 10)
	arm=$(expr $arm + $r \* $r \* $r)
	n=$(expr $n / 10)
done
echo $arm
if [ $arm -eq $temp ]; then
	echo "Armstrong"
else
	echo "Not Armstrong"
fi
