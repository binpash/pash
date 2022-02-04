DATA=final.csv
rm -f $DATA
replace_string() {
    sed -i -e 's/'$1'/'"$2"'/g' .tmp
}

prepare_run_data() {
    DATA_FILE=$1.data.csv
    FILES=$1.log.dat
    echo $1.tmp
    rm -f $DATA_FILE $FILES

    # gather all results
    find $1 -name "*.res" -type f > $FILES

    # read each line of the file
    while read p; do
        lines=$(cat $p | wc -l)
        if [[ $lines -gt 3 ]]; then
            file="$(tail -n +3 $p)"
        else
            file="$(tail -n +2 $p)"
        fi
        echo "$file" > .tmp
        bench=$(echo $p | awk -F '/' '{print $3}')
        mode=$(echo $p | awk -F'/' '{print $2}')
        # read the contents of each execution file
        while read l; do
            #l=$(echo $l | sed 's/ //g')
            echo $l | grep --quiet :
            res=$?
            if [[ $res == 1 ]]; then
                perf=$(echo $l | grep -Eo '[0-9]+.[0-9]+$')
                script=$(echo $l | sed -e 's/'$perf'//g')
            else
                # get script name and performance
                # strip the .sh and get fetch the script name
                script=$(echo $l | awk -F ':' '{print $1}' | sed 's/...$//')
                # get the execution time
                perf=$(echo $l | awk -F ':' '{print $2}')
            fi
            echo $bench,$script,$mode,$perf | sed 's/ //g'>> $DATA_FILE
        done < .tmp
    done < $FILES
    sort $DATA_FILE > $1.tmp
    rm $DATA_FILE rm -f $FILES
}

cd eval_results
prepare_run_data run1
prepare_run_data run2
prepare_run_data run3
# merge all the results
paste run1.tmp run2.tmp run3.tmp -d , | sed -s  's/,/ /g' | awk '{print $1,$2,$3, $4+$8+$12}' | awk ' {print $1','$2','$3','$4/3}' | tr ' ' ',' > .tmp
# cleanup
replace_string dependency_untangling for-loops 
replace_string nlp NLP
replace_string oneliners Classics
replace_string unix50 Unix50
replace_string analytics-mts COVID-mts
replace_string web-index WebIndex
replace_string max-temp AvgTemp
replace_string temp-analytics AvgTemp
replace_string Genomics_Computation Genomics
replace_string Program_Inference ProgInf
replace_string blish_no_prof_no_du 'pash_jit -prof -par_pipe'
replace_string blish_no_prof 'pash_jit -prof'
replace_string pash pash_aot
replace_string blish pash_jit
perf=''
# calculate the ratios
while read p; do
    # is this the bash entry
    echo $p | grep --quiet bash
    res=$?
    # fetch the performance
    if [[ $res == 0 ]]; then
        perf=$(echo $p | awk -F ',' '{print $4}')
    fi
    # get the bench
    bench=$(echo $p | awk -F ',' '{print $1}')
    # get the script
    script=$(echo $p | awk -F ',' '{print $2}')
    # get the mode
    mode=$(echo $p | awk -F ',' '{print $3}')
    # get the time of the pash/blish configs
    current_perf=$(echo $p | awk -F ',' '{print $4}')
    # calculate the ratio
    if [[ $res == 0 ]]; then
        ratio=$perf
    else
        ratio=$(echo "$perf $current_perf" | awk '{print $1/$2}' )
    fi
    # replace the pash/blish time with the ratio
    echo $bench,$script,$mode,$ratio >> $DATA
done < .tmp
rm -f .tmp
mv $DATA ..
