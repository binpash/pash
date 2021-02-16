rm -f $PASH_TOP/evaluation/scripts/input/error.txt
rm -f $PASH_TOP/evaluation/scripts/input/succ.txt
# strip the first column
mkdir -p $PASH_TOP/evaluation/scripts/input/generated
cut -f1 -d" " --complement  $PASH_TOP/evaluation/scripts/input/likely-longest-pipelines.txt > $PASH_TOP/evaluation/scripts/input/generated/generated.file

cd $PASH_TOP/evaluation/scripts/input/generated
# we could read the file iteratively with IFS, but the environment was affected
# idk maybe Im bad
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
