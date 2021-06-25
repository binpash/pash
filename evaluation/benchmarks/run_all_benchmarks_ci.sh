#!/bin/bash
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
## This script is necessary to ensure that sourcing happens with bash
source run.seq.sh
source run.par.sh
# total tests
total=0
# number of tests that passed
passed=0
compare_outputs(){
    dir=$1
    outputs=$(ls $dir | grep "seq" | sed 's/.seq.out$//')
    for out in $outputs;
    do
        seq_output="${dir}/${out}.seq.out"
        pash_output="${dir}/${out}.par.out"
        res=$(diff -q "$seq_output" "$pash_output")
        if [[ "${res}" -eq "" ]]; then
            passed=$((passed + 1))
        fi
        total=$((total + 1))
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
b=$($1      | grep -v executing | sed 's\.sh:\bash\g' |  sed 's\,\.\g' | awk '{ print $1 "," $2}')
labels="group,Bash"
for conf in "${configurations[@]}"; do
    for n_in in "${n_inputs[@]}"; do
        # on each run, clean all the res files 
        rm -f $1/*.res
        # re-export the new config
        export PASH_FLAGS="${conf} -w ${n_in}"
        # append the new labels for the plot
        labels="${labels},${conf}_${n_in}"
        # execute the pash with the new config
        p=$($1_pash | grep -v executing | sed 's\.sh:\pash\g' |  sed 's\,\.\g' | awk '{ print $1 "," $2}' | awk '{sub(/[^,]*/,"");sub(/,/,"")} 1')
        # store the results
        EXEC+=("$p")
    done
done
# concat all the results and merge them to create the final data for plotting
labels=$(echo $labels | sed 's\--\\g' | sed -e 's/ /_/g')
res="$b"
for i in "${EXEC[@]}"
do
    res=$(paste -d'@' <(echo "$res") <(echo "$i"))
done
# write the labels to the file
echo "$labels" > results.time
# write the data formatted
echo -e "$res" | sed 's\@\,\g' >> results.time
# compare the results
compare_outputs "$1/outputs"
# this is going to be written on the UI output log
cat results.time
# this is going to be written on the UI output log / CLI output
echo "Summary: ${passed}/${total} tests passed."
