#!/bin/bash

# exit when any command fails
#set -e

# Download inputs
# 						~ Use two sources: University (fast) 
# 					             			Long term storage (slow)
# 						~ And two sizes:  Small & quick
# 	             Full size

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs
hdfs dfs -mkdir /file-enc
                                                  

# download the initial pcaps to populate the whole dataset
if [ ! -d ${IN}/pcap_data ]; then
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

  # generates small inputs
  mkdir -p pcap_data_small/
  PCAP_DATA_FILES=1
  for (( i = 1; i <= $PCAP_DATA_FILES; i++ )) do
      for j in pcaps/*;do
          n=$(basename $j)
          cat $j > pcap_data_small/pcap${i}_${n}; 
      done
  done
  echo "Pcaps_small Generated"
fi
