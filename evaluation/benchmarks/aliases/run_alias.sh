# parses the generated.file, and creates a log of the commands that were executed
# successfully (succ.txt) and the failed ones (err.txt)

cd $PASH_TOP/evaluation/scripts/input/
# we could read the file iteratively with IFS, but the environment was affected
IFS=$'\r\n' GLOBIGNORE='*' command eval  'cmd_array=($(cat generated.file))'
lc=$(cat generated.file | wc -l)
for i in $(seq 0 $lc)
do
    # get the entry from the array
    p=${cmd_array[$i]}
    # add a timeout to our script
    timeout --signal=SIGINT 50s /bin/bash -e $p >> /dev/null 2>&1  #./cmd.sh #eval "bash ./cmd.sh"
    ## get status ##
    status=$?
    if [ $status -eq 0 ]; then
        echo $p >> $PASH_TOP/evaluation/scripts/input/succ.txt
    else
        echo $p >> $PASH_TOP/evaluation/scripts/input/err.txt
    fi
    if ! ((i % 100)); then
        echo $i
    fi
done
echo "Done"
