#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

generate_jpg_manifest() {
    local manifest_path="../jpg_list.txt"
    if [ -d "jpg_full/jpg" ]; then
        find "jpg_full/jpg" -maxdepth 1 -type f -printf '%f\n' | sort > "$manifest_path"
    fi
}

# download the input for the nginx logs and populate the dataset
if [ ! -d log_data_heavy ]; then
    wget https://atlas-group.cs.brown.edu/data/nginx.zip
    unzip nginx.zip 
    rm nginx.zip
    # generating full analysis logs
    mkdir -p log_data_heavy
    LOG_DATA_FILES=1
    for (( i = 1; i <=$LOG_DATA_FILES; i++)) do
        for j in nginx-logs/*; do
            n=$(basename $j)
            # cal the times need to generate 500M data
            iters=$((500 * 1024 * 1024 / $(stat -c%s "$j")))
            for (( k = 1; k <= $iters; k++ )) do
                cat $j >> log_data_heavy/log${i}0_500M_${n}_500M.log.log;
            done
        done
    done
    echo "Nginx logs Generated"

	# # generating small analysis logs
    # mkdir -p log_data_small
	# LOG_DATA_FILES=6
    # for (( i = 1; i <=$LOG_DATA_FILES; i++)) do
    #     for j in nginx-logs/*;do
    #         n=$(basename $j)
    #         cat $j > log_data_small/log${i}_${n}.log; 
    #     done
    # done
    # echo "Nginx logs (small) Generated"
fi


if [ ! -d pcap_data_heavy ]; then
  wget https://atlas-group.cs.brown.edu/data/pcaps.zip
  unzip pcaps.zip
  rm pcaps.zip
  # generates 20G
  mkdir -p pcap_data_heavy/
  PCAP_DATA_FILES=64
    for j in pcaps/*;do
        n=$(basename $j)
        iters=$((20 * 1024 * 1024 * 1024 / $(stat -c%s "$j")))
        for (( k = 1; k <= $iters; k++ )) do
            cat $j >> pcap_data_heavy/500M_${n};
        done
    done
  echo "Pcaps Generated"
fi

if [ ! -d jpg_full/jpg ]; then
  wget https://atlas-group.cs.brown.edu/data/full/jpg.zip -O jpg_full.zip
  unzip jpg_full.zip -d jpg_full
  rm jpg_full.zip
  echo "JPG Generated"
fi

generate_jpg_manifest

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
