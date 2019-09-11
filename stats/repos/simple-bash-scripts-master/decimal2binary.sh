#!/bin/bash

for ((i = 32; i >= 0; i--)); do
	r=$((2 ** $i))
	Probablity+=($r)
done

[[ $# -eq 0 ]] && {
	echo -e "Usage \n \t $0 numbers"
	exit 1
}

echo -en "Decimal\t\tBinary\n"
for input_int in $@; do
	s=0
	test ${#input_int} -gt 11 && {
		echo "Support Upto 10 Digit number :: skiping \"$input_int\""
		continue
	}

	printf "%-10s\t" "$input_int"

	for n in ${Probablity[@]}; do

		if [[ $input_int -lt ${n} ]]; then
			[[ $s == 1 ]] && printf "%d" 0
		else
			printf "%d" 1
			s=1
			input_int=$(($input_int - ${n}))
		fi
	done
	echo -e
done
