#!/bin/bash
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
## This script is necessary to ensure that sourcing happens with bash
source run.seq.sh
source run.par.sh

compare_outputs(){
  dir=$1
  outputs=$(ls $dir | grep "seq" | sed 's/.seq.out$//')
  for out in $outputs;
  do
    seq_output="${dir}/${out}.seq.out"
    pash_output="${dir}/${out}.par.out"
    diff -q "$seq_output" "$pash_output"
  done
}

#oneliners
#oneliners_pash
#
#compare_outputs "oneliners/outputs"
#
#unix50
#unix50_pash
#
#compare_outputs "unix50/outputs"
#
#poets
#poets_pash
#
#compare_outputs "poets/outputs"
#
#web-index
#web-index_pash
#
#compare_outputs "web-index/outputs"
EXPERIMENTAL=1
if [ "$EXPERIMENTAL" -eq 1 ]; then
    configurations=(
        # "" # Commenting this out since the tests take a lot of time to finish
        "--r_split"
        "--dgsh_tee"
        "--r_split --dgsh_tee"
        # "--speculation quick_abort"
    )
else
    configurations=(
        ""
    )
fi


n_inputs=(
    2
    8
    16
)

EXEC=()

# cleanup
rm -f $1/*.res
# run bash
b=$($1      | grep -v executing | sed 's\:\\g'|  sed 's\,\.\g' | awk '{ print $1 "," $2}')

labels="date,hash,conf,script,time,id"
hash=$(git rev-parse --short HEAD)
ts=$(date +"%Y-%m") 
b=$(echo "$b" | sed -e "s/^/$ts,$hash,/")
b=$(echo "$b" | sed -e "s/$/,bash/" | sed 's\--\\g' | sed 's\ \\g')

echo "$labels" > results.time
echo "$b" >> results.time




rm -rf output_res
for conf in "${configurations[@]}"; do
    for n_in in "${n_inputs[@]}"; do
        rm -f $1/*.res
        #echo "|-- Executing with pash --width ${n_in} ${conf}..."
        export PASH_FLAGS="${conf} -w ${n_in}"
        p=$($1_pash | grep -v executing | sed 's\:\\g' | sed 's\,\.\g' | awk '{ print $1 "," $2}' )
        
        p=$(echo "$p" | sed -e "s/^/$ts,$hash,/")
        p=$(echo "$p" | sed -e "s/$/,${conf}_${n_in}/" | sed 's\--\\g' | sed 's\ \\g')
        
        echo "$p"
        EXEC+=("$p")
    done
done
labels=$(echo $labels | sed 's\--\\g' | sed -e 's/ /_/g')
for i in "${EXEC[@]}"
do
    echo "$i" >> results.time
done
compare_outputs "$1/outputs"
cat results.time
