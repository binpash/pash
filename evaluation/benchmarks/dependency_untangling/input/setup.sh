#!/bin/bash

# exit when any command fails
#set -e

IN=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/
OUT=$PASH_TOP/evaluation/benchmarks/dependency_untangling/output/
IN_NAME=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/100G.txt

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
    rm -rf ${IN}/bio
    rm -rf ${IN}/output
    rm -rf ${OUT}
    exit 
fi

setup_dataset() {
  if [ "$1" == "--small" ]; then
      LOG_DATA_FILES=6
      WAV_DATA_FILES=20
      NODE_MODULE_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/node_modules.zip
      BIO_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/bio.zip
      JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/jpg.zip
      PCAP_DATA_FILES=1
  else
      LOG_DATA_FILES=84
      WAV_DATA_FILES=120
      NODE_MODULE_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/node_modules.zip
      BIO_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/bio.zip
      JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/jpg.zip
      PCAP_DATA_FILES=15
  fi
  
  if [ ! -d ${IN}/wav ]; then
      wget http://pac-n4.csail.mit.edu:81/pash_data/wav.zip
      unzip wav.zip && cd wav/
      for f in *.wav; do
          FILE=$(basename "$f")
          for (( i = 0; i <= $WAV_DATA_FILES; i++)) do
              echo copying to $f$i.wav
              cp $f $f$i.wav
          done
      done
      echo "WAV Generated"
  fi
  
  if [ ! -d ${IN}/jpg ]; then
      cd ${IN}
      wget $JPG_DATA_LINK
      unzip jpg.zip
      echo "JPG Generated"
      rm -rf ${IN}/jpg.zip
  fi
  
  # download the input for the nginx logs and populate the dataset
  if [ ! -d ${IN}/log_data ]; then
      cd $IN
      wget http://pac-n4.csail.mit.edu:81/pash_data/nginx.zip
      unzip nginx.zip 
      rm nginx.zip
      # generating analysis logs
      mkdir -p ${IN}/log_data
      for (( i = 1; i <=$LOG_DATA_FILES; i++)) do
          for j in nginx-logs/*;do
              n=$(basename $j)
              cat $j > log_data/log${i}_${n}.log; 
          done
      done
      echo "Logs Generated"
  fi
  
  if [ ! -d ${IN}/bio ]; then                                                  
      if [ "$1" = "--small" ]; then
          # download the Genome loc file
          wget $BIO_DATA_LINK 
          unzip bio.zip
          cd bio
          wget http://pac-n4.csail.mit.edu:81/pash_data/Gene_locs.txt
          wget http://pac-n4.csail.mit.edu:81/pash_data/small/100G.txt
          cd ..
          rm bio.zip
      else
          mkdir ${IN}/bio                                                          
          cd ${IN}/bio                                                             
          # download the file containing the links for the dataset                 
          wget http://pac-n4.csail.mit.edu:81/pash_data/100G.txt                   
          # download the Genome loc file                                           
          wget http://pac-n4.csail.mit.edu:81/pash_data/Gene_locs.txt              
          # start downloading the real dataset  
          IN_NAME=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/bio/100G.txt
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
      fi                                                                           
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
      for (( i = 1; i <= $PCAP_DATA_FILES; i++ )) do
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
      wget $NODE_MODULE_LINK
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
      if [ "$1" = "--small" ]; then
          head -n 20 packages > p
          mv p  packages
      fi
      echo "Package datset downloaded"
  fi
}

source_var() {
  export IN=
}
