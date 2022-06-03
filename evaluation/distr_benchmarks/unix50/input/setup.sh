#!/bin/bash

#set -e

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

## FIXME: These inputs are already 1G when downloaded
## FIXME: Also, wget is not silent like curl in the other setup scripts.

inputs=(
1 10 11 12 2 3 4 5 6 7 8 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9
)

if [[ "$1" == "-c" ]]; then
    for input in ${inputs[@]}
    do
        rm -f "${input}.txt"
    done
    exit
fi

setup_dataset() {
    # Put files in hdfs
    hdfs dfs -mkdir /unix50
    
    # generate small inputs 
    # if [ "$#" -eq 1 ] && [ "$1" = "--small" ]; then
    #   if [ ! -d ./small ]; then                                                          
    #     echo "Generating small-size inputs"                                             
    #     # FIXME PR: Do we need all of them?                                             
    #     curl -sf 'http://pac-n4.csail.mit.edu:81/pash_data/small/unix50.zip' > unix50.zip
    #     unzip unix50.zip                                                                 
    #     rm -f unix50.zip                                                                 
    #   fi
    #   hdfs dfs -put small /unix50/small                                                                              
    #   return 0
    # fi
  
    for input in ${inputs[@]}
    do
        if [ ! -f "${input}.txt" ]; then
            wget "http://ndr.md/data/unix50/${input}.txt"
            "$PASH_TOP/scripts/append_nl_if_not.sh" "${input}.txt"
        fi
        hdfs dfs -put $file /unix50/$file
    done

    # increase the original input size 10x
    if [ "$#" -eq 1 ] && [ "$1" = "--extended" ]; then
        EXTENDED_INPUT_DIR="extended_input/"
        mkdir -p $EXTENDED_INPUT_DIR
        for file in *.txt; do
            rm $EXTENDED_INPUT_DIR/$file
            for (( i = 0; i < 10; i++ )); do
                cat $file >> $EXTENDED_INPUT_DIR/temp.txt
            done
        done
        hdfs dfs -put $EXTENDED_INPUT_DIR /unix50/$EXTENDED_INPUT_DIR
        rm -rf $EXTENDED_INPUT_DIR
    fi
}

source_var() {
  if [[ "$1" == "--extended" ]]; then
    export IN_PRE=/unix50/extended_input
  else
    export IN_PRE=/unix50
  fi
}
