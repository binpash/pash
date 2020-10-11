#!/bin/bash
echo "Enter a number :"
read Binary
if [ $Binary -eq 0 ]; then
	echo "Enter a valid number "
	return
else
	while [ $Binary -ne 0 ]; do
		Bnumber=$Binary
		Decimal=0
		power=1
		while [ $Binary -ne 0 ]; do
			rem=$(expr $Binary % 10)
			Decimal=$((Decimal + (rem * power)))
			power=$((power * 2))
			Binary=$(expr $Binary / 10)
		done
		echo " $Decimal"
	done
fi
