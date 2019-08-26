#!/bin/bash
echo "Enter The Number"
read a
fact=1
while [ $a -ne 0 ]; do
	fact=$(expr $fact \* $a)
	a=$(expr $a - 1)
done
echo $fact
