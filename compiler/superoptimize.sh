#!/bin/bash
DIGITS=3
TIMEFORMAT="RadhasTimer %${DIGITS}R"

echo "Superotmizer run!"

NPROC=$(nproc) # is this max width? or pick some after that
echo $NPROC 

COUNTER=1
widths=() # determine what intervals 
for ((i=1; i<=$NPROC; i++)); do
   widths[i]=$COUNTER
   COUNTER=$(($COUNTER+1))
done

echo "${widths[@]}"

intensive_commands=("sort" "grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'")
input=("1M.txt" "10M.txt" "100M.txt" "1G.txt")

echo "Begin analysis"
for f in "${intensive_commands[@]}"; do
    echo $f
    for inp in "${input[@]}"; do 
        echo $inp
        IN=../evaluation/scripts/input/${inp}
        echo "cat $IN | $f" > test-script.sh
        for width in "${widths[@]}"; do 
            echo $width > percommand.sh
            (time $PASH_TOP/pa.sh -w $width -d 0 test-script.sh > out.txt) 2>&1 | grep RadhasTimer | sed 's/RadhasTimer //' > ./result.sh
            RESULTS="$(cat result.sh)w${width}"
            echo $RESULTS >> ./results.sh
        done 
        echo "${f}, ${inp}: $(cat ./results.sh | sort -k1 -tw -n -r | head -1)" >> optimal-config.sh
        rm results.sh
    done
done

# parse script, find most intensive command -> look into optimal-config doc to find width 
# add space component; have to determine though what is "max" space (use another flag?) as well as how much space each eager file takes up
        