set -o pipefail
# parses the generated.file, and creates a log of the commands that were executed
# successfully (succ.txt) and the failed ones (err.txt)
cd $PASH_TOP/evaluation/scripts/input/
rm -f succ.txt
rm -f err.txt
# we could read the file iteratively with IFS, but the environment was affected
IFS=$'\r\n' GLOBIGNORE='*' command eval  'cmd_array=($(cat generated.file))'
lc=$(cat generated.file | wc -l)
for i in $(seq 0 $lc)
do
    # get the entry from the array
    p=${cmd_array[$i]}
    # remove start/end quotes
    v=$(echo $p  | sed -e 's/^"//' -e 's/"$//')
    # write command to file
    echo "set -o pipefail" > .tmp.sh
    echo "$v > /dev/null">> .tmp.sh
    # add a timeout to our script
    timeout -k 50 --signal=SIGINT 50s /bin/bash -e ./.tmp.sh > /dev/null 2>&1
    if [ $? -ne 0 ]; 
    then
        echo $v >> err.txt
    else
        echo $v >> succ.txt
    fi
    if ! ((i % 100)); then
        echo $i
    fi



done
echo "Done"
