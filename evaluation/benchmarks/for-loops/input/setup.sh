#!/bin/bash

# exit when any command fails
#set -e

IN=$PASH_TOP/evaluation/benchmarks/for-loops/input/
OUT=$PASH_TOP/evaluation/benchmarks/for-loops/output/
IN_NAME=$PASH_TOP/evaluation/benchmarks/for-loops/input/bio/100G.txt
if [ "$1" == "-c" ]; then
	rm -rf ${IN}/jpg
	rm -rf ${IN}/log_data
	rm -rf ${IN}/wav
	rm -rf ${IN}/nginx-logs
	rm -rf ${IN}/node_modules
	rm -rf ${IN}/pcap_data
    rm -rf ${IN}/pcaps
    rm -rf ${IN}/packages
    rm -rf ${IN}/mir-sa
    rm -rf ${IN}/deps
	rm -rf ${OUT}
	exit 
fi

if [ ! -d ${IN}/wav ]; then
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

# download the input for the nginx logs and populate the dataset
if [ ! -d ${IN}/log_data ]; then
	cd $IN
	wget http://pac-n4.csail.mit.edu:81/pash_data/nginx.zip
	unzip nginx.zip 
    rm nginx.zip
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

if [ ! -d ${IN}/bio ]; then
	mkdir ${IN}/bio
	cd ${IN}/bio
	# download the file containing the links for the dataset
	wget http://pac-n4.csail.mit.edu:81/pash_data/100G.txt
	# download the Genome loc file
	wget http://pac-n4.csail.mit.edu:81/pash_data/Gene_locs.txt
	# start downloading the real dataset
	cat ${IN_NAME} |while read s_line;
	do
		echo ${IN_NAME}
		sample=$(echo $s_line |cut -d " " -f 2);
		if [[ ! -f $sample ]]; then
			pop=$(echo $s_line |cut -f 1 -d " ");
			link=$(echo $s_line |cut -f 3 -d " ");
			wget -O "$sample".bam  "$link"; ##this part can be adjusted maybe
		fi
	done;
    echo "Genome data downloaded"
fi

# download the initial pcaps to populate the whole dataset
if [ ! -d ${IN}/pcap_data ]; then
	cd $IN
	wget http://pac-n4.csail.mit.edu:81/pash_data/pcaps.zip
	unzip pcaps.zip
    rm pcaps.zip
	mkdir ${IN}/pcap_data/
	# generates 20G
	for i in {1..15};do
		for j in ${IN}/pcaps/*;do
			n=$(basename $j)
			cat $j > pcap_data/pcap${i}_${n}; 
		done
	done
	echo "Pcaps Generated"
fi 

# download the modules for the Mir static analyses
if [ ! -d ${IN}/node_modules ]; then
	cd $IN
	wget http://pac-n4.csail.mit.edu:81/pash_data/node_modules.zip
	unzip node_modules.zip 
    rm node_modules.zip
	# download the specific mir version
	wget http://pac-n4.csail.mit.edu:81/pash_data/mir-sa.zip
	unzip mir-sa.zip
    rm mir-sa.zip
	echo "Node modules generated"
fi

# download the packages for the package building
if [ ! -f ${IN}/packages ]; then
	cd $IN
	wget http://pac-n4.csail.mit.edu:81/pash_data/packages
	echo "Package datset downloaded"
fi
