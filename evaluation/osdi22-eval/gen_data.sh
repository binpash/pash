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
        # echo "Fixing: $p"
        python3 ../prep_temp.py $p > .tmp
        # break        
        # lines=$(cat $p | wc -l)
        # if [[ $lines -gt 3 ]]; then
        #     echo "Head1:"
        #     tail -n +2 $p | head
        #     file="$(tail -n +2 $p)"
        # else
        #     echo "Head2:"
        #     tail -n +2 $p | head
        #     file="$(tail -n +2 $p)"
        # fi
        # file="$(tail -n +2 $p)"
        # echo "file: $file"
        # echo "$file" > .tmp
        bench=$(echo $p | awk -F '/' '{print $3}')
        mode=$(echo $p | awk -F'/' '{print $2}')
        # echo "Bench: $bench, mode: $mode"
        # if [[ $bench == max-temp ]]; then
        # cat .tmp | sed -E 's/^([a-zA-Z_0-9\-]+):.*([0-9]+.[0-9]+\n)$/\1\t\2/g' #| cut -f 1 
        # fi
        # read the contents of each execution file
        while read l; do
            #l=$(echo $l | sed 's/ //g')
            echo $l | grep --quiet :
            res=$?
            if [[ $res == 1 ]]; then
                perf=$(echo $l | grep -Eo '[0-9]+.[0-9]+$')
                # echo "Perf 1: $perf"
                script=$(echo $l | sed -e 's/'$perf'//g')
            else
                # get script name and performance
                # strip the .sh and get fetch the script name
                script=$(echo $l | awk -F ':' '{print $1}' | sed 's/...$//')
                # get the execution time
                perf=$(echo $l | awk -F ':' '{print $2}')
                # echo "Perf 2: $perf"
            fi
            echo $bench,$script,$mode,$perf | sed 's/ //g'>> $DATA_FILE
        done < .tmp
    done < $FILES
    sort $DATA_FILE > $1.tmp
    rm $DATA_FILE rm -f $FILES
}

cd eval_results
prepare_run_data run
# merge all the results
cat run.tmp | sed -s  's/,/ /g' | awk '{print $1,$2,$3,$4}' | awk ' {print $1','$2','$3','$4}' | tr ' ' ',' > .tmp
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
replace_string pash_jit_no_prof_no_du 'pash_jit -prof -par_pipe'
replace_string pash_jit_no_prof 'pash_jit -prof'
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

# in docker container, we are running with the CI
if [ -f /.dockerenv ]; then
    exit 0
fi
cd ..
# replace all the lines that are not needed in figure5
sed 's/for-loops,AurPkg,pash_aot,.*/for-loops,AurPkg,pash_aot,0/g' $DATA > data_final.csv
sed -i 's/for-loops,FileEnc1,pash_aot,.*/for-loops,FileEnc1,pash_aot,0/g' data_final.csv 
sed -i 's/for-loops,FileEnc2,pash_aot,.*/for-loops,FileEnc2,pash_aot,0/g' data_final.csv 
sed -i 's/for-loops,LogAnalysis1,pash_aot,.*/for-loops,LogAnalysis1,pash_aot,0/g' data_final.csv
sed -i 's/for-loops,LogAnalysis2,pash_aot,.*/for-loops,LogAnalysis2,pash_aot,0/g' data_final.csv
sed -i 's/for-loops,MediaConv1,pash_aot,.*/for-loops,MediaConv1,pash_aot,0/g' data_final.csv
sed -i 's/for-loops,MediaConv2,pash_aot,.*/for-loops,MediaConv2,pash_aot,0/g' data_final.csv
sed -i 's/for-loops,ProgInf,pash_aot,.*/for-loops,ProgInf,pash_aot,0/g' data_final.csv
sed -i 's/for-loops,Genomics,pash_aot,.*/for-loops,Genomics,pash_aot,0/g' data_final.csv
