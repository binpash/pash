#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

# download the input for the nginx logs and populate the dataset
if [ ! -d log_data ]; then
    wget https://atlas-group.cs.brown.edu/data/nginx.zip
    unzip nginx.zip 
    rm nginx.zip
    # generating full analysis logs
    mkdir -p log_data
	LOG_DATA_FILES=84
    for (( i = 1; i <=$LOG_DATA_FILES; i++)) do
        for j in nginx-logs/*;do
            n=$(basename $j)
            cat $j > log_data/log${i}_${n}.log; 
        done
    done
    echo "Nginx logs Generated"

	# generating small analysis logs
    mkdir -p log_data_small
	LOG_DATA_FILES=6
    for (( i = 1; i <=$LOG_DATA_FILES; i++)) do
        for j in nginx-logs/*;do
            n=$(basename $j)
            cat $j > log_data_small/log${i}_${n}.log; 
        done
    done
    echo "Nginx logs (small) Generated"
fi


if [ ! -d pcap_data ]; then
  wget https://atlas-group.cs.brown.edu/data/pcaps.zip
  unzip pcaps.zip
  rm pcaps.zip
  # generates 20G
  mkdir -p pcap_data/
  PCAP_DATA_FILES=15
  for (( i = 1; i <= $PCAP_DATA_FILES; i++ )) do
      for j in pcaps/*;do
          n=$(basename $j)
          cat $j > pcap_data/pcap${i}_${n};
      done
  done
  echo "Pcaps Generated"
fi

#   # generates small inputs
#   mkdir -p pcap_data_small/
#   PCAP_DATA_FILES=1
#   for (( i = 1; i <= $PCAP_DATA_FILES; i++ )) do
#       for j in pcaps/*;do
#           n=$(basename $j)
#           cat $j > pcap_data_small/pcap${i}_${n}; 
#       done
#   done
#   echo "Pcaps_small Generated"
# fi
