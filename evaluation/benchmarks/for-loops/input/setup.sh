#!/bin/bash

# exit when any command fails
#set -e


IN=$PASH_TOP/evaluation/benchmarks/for-loops/input/
OUT=$PASH_TOP/evaluation/benchmarks/for-loops/output/
IN_NAME=${IN_N:-PASH_TOP/evaluation/benchmarks/foor-loops/input/100G.txt}
if [[ "$1" == "-c" ]]; then
    rm -rf ${IN}/jpg
    rm -rf ${IN}/rtf
    rm -rf ${IN}/wav
    rm -rf ${IN}/linux
    rm -rf ${IN}/logs
    rm -rf ${IN}/bio
    rm -rf ${OUT}
    exit 
fi

if [[ ! -d ${IN}/wav  ]]; then
    mkdir -p ${IN}/wav
    cd ${IN}/wav
    wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_1MG.wav
    wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_2MG.wav
    wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_5MG.wav
    wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_10MG.wav

    for f in *.wav; do
        FILE=$(basename "$f")
        for i in {1..120}; do 
            echo copying to $f$i.wav
            cp $f $f$i.wav
        done
    done
    echo "WAV Generated"
fi
if [ ! -d ${IN}/jpg ]; then
    mkdir -p ${IN}/jpg
    cd ${IN}/jpg
    curl ftp://ftp.inrialpes.fr/pub/lear/douze/data/jpg1.tar.gz | tar -xz
    mv jpg/* . && rm -r jpg/

    mkdir -p tmpd
    for filename in *.jpg; do                                                     
        cp $filename tmpd/${filename}_copy.jpg   
    done

    for filename in tmpd/*.jpg; do                                                     
        cp $filename ${IN}/jpg
    done                                                                          
    echo "JPGs copied"
    rm -rf ${IN}/jpg/tmpd
fi

if [ ! -d ${IN}/log_data ]; then
    cd $IN
    # generating analysis logs
    mkdir -p ${IN}/log_data
    for i in {1..84};do 
        for j in nginx-logs/*;do
            n=$(basename $j)
            cat $j > log_data/log${i}_${n}.log; 
        done
    done
    echo "Logs Generated"
fi
#cat ${IN}/${IN_NAME} |while read s_line;
#	do
#    
#    sample=$(echo $s_line |cut -d " " -f 2);
#    if [[ ! -f $sample ]]; then
#        pop=$(echo $s_line |cut -f 1 -d " ");
#        link=$(echo $s_line |cut -f 3 -d " ");
#        wget -O "$PW/$sample".bam  "$link"; ##this part can be adjusted maybe
#    fi
#done;
#

if [ ! -d ${IN}/pcap_data ]; then
    mkdir ${IN}/pcap_data/
    cd ${IN}/pcap_data/
    # generates 20G
    for i in {1..15};do
        for j in ${IN}/pcaps/*;do
            n=$(basename $j)
            cat $j > pcap_data/pcap${i}_${n}; 
        done
    done
    echo "Pcaps Generated"
fi 

    install() {

        file=$1
        while IFS= read -r package
        do
            echo "Installing $package..."
            npm list $package >/dev/null 2>/dev/null || npm install $package >/dev/null 2>error.log
            #if grep -q "ERR!" error.log; then
            #  echo "$1" >> ../$package.results/errOnInstall 
            #fi
            #rm error.log
        done < "$file"
    }

if [ ! -d ${IN}/node_modules ]; then
    npm i -g @andromeda/mir-sa
    mv node_modules mir-sa
    install ${IN}/mir-packages.txt
    echo "Node modules downloaded"
fi

if [ ! -f ${IN}/packages ]; then
    cd ${IN}
    wget https://aur.archlinux.org/packages.gz
    gunzip packages.gz
    head -n 150 packages > t
    mv t packages;
    echo "Packages Extracted"
fi
