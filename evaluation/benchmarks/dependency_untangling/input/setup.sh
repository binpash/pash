#!/bin/bash

# exit when any command fails
#set -e

IN=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/
OUT=$PASH_TOP/evaluation/benchmarks/dependency_untangling/output/
IN_NAME=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/100G.txt

## Import the build library
. "$IN/build_lib.sh"

if [ "$1" == "-c" ]; then
    rm -rf "${IN}/jpg"
    rm -rf "${IN}/log_data"
    rm -rf "${IN}/wav"
    rm -rf "${IN}/nginx-logs"
    rm -rf "${IN}/node_modules"
    rm -rf "${IN}/pcap_data"
    rm -rf "${IN}/pcaps"
    rm -rf "${IN}/packages"
    rm -rf "${IN}/mir-sa"
    rm -rf "${IN}/deps"
    rm -rf "${IN}/bio"
    rm -rf "${IN}/output"
    rm -rf "${OUT}"
    exit 
fi



## Q: Can these checks be generated automatically? This would be great if
##    the user just ran the command, and then if it succeeded, the test is generated.
wav_step_1_done_check()
{
    local prefix="wav/file_example_WAV"
    files_exist_done_check "${prefix}_1MG.wav.kernel" "${prefix}_2MG.wav.kernel" "${prefix}_5MG.wav.kernel" "${prefix}_10MG.wav.kernel"
    return $?
}

## Q: Can we automatically check that this is idempotent?
## For example, in the step below there were 2 non-idempotence issues:
## - wget downloads wav.zip.2 if wav.zip already exists, so we need to use -O flag
## - wav files need to be saved with .kernel suffix to make step 2 idempotent
wav_step_1()
{
    curl -C - -o wav.zip http://pac-n4.csail.mit.edu:81/pash_data/wav.zip
    unzip wav.zip
    local prefix="wav/file_example_WAV"
    ## Necessary so that the iteration in step 2 is idempotent
    for f in "${prefix}_1MG.wav" "${prefix}_2MG.wav" "${prefix}_5MG.wav" "${prefix}_10MG.wav"; do
        mv $f $f.kernel
    done
}
export -f wav_step_1_done_check
export -f wav_step_1

wav_step_2_done_check()
{
    local prefix="wav/file_example_WAV"
    for i in $(seq 0 "$WAV_DATA_FILES"); do
        if ! files_exist_done_check "${prefix}_1MG.wav$i.wav" "${prefix}_2MG.wav$i.wav" "${prefix}_5MG.wav$i.wav" "${prefix}_10MG.wav$i.wav"; then
            return 1
        fi
    done
    echo "Done"
    return 0
}

wav_step_2()
{
    for f in wav/*.kernel; do
        for (( i = 0; i <= $WAV_DATA_FILES; i++)) do
            echo copying to "$base_f$i.wav"
            base_f=wav/$(basename "$f" .kernel)
            cp "$f" "$base_f$i.wav"
        done
    done
}
export -f wav_step_2_done_check
export -f wav_step_2


jpg_step()
{
    curl -C - -o jpg.zip $JPG_DATA_LINK
    unzip jpg.zip
    rm -rf ${IN}/jpg.zip
}

jpg_step_done_check()
{
    number_of_files_in_dir $JPG_NUMBER jpg
}
export -f jpg_step
export -f jpg_step_done_check

nginx_logs_step_1()
{
    curl -C - -o nginx.zip http://pac-n4.csail.mit.edu:81/pash_data/nginx.zip
    unzip nginx.zip
    rm nginx.zip
}

nginx_logs_step_1_done_check()
{
    local prefix="nginx-logs/log"
    for i in $(seq 0 7); do
        if ! files_exist_done_check "${prefix}$i"; then
            return 1
        fi
    done
    return $?
}

export -f nginx_logs_step_1
export -f nginx_logs_step_1_done_check

nginx_logs_step_2()
{
    # generating analysis logs
    mkdir -p ${IN}/log_data
    for (( i = 1; i <=$LOG_DATA_FILES; i++)) do
        for j in nginx-logs/*;do
            n=$(basename $j)
            cp $j  log_data/log${i}_${n}.log; 
        done
    done
}


nginx_logs_step_2_done_check()
{
    local prefix="log_data/log"
    for j in $(seq 0 7); do
        for i in $(seq 1 "$LOG_DATA_FILES"); do
            if ! files_exist_done_check "${prefix}${i}_log${j}.log"; then
                return 1
            fi
        done
    done
    echo "Done"
    return 0
}

export -f nginx_logs_step_2_done_check
export -f nginx_logs_step_2


setup_dataset() {
  if [ "$1" == "--small" ]; then
      export LOG_DATA_FILES=6
      export WAV_DATA_FILES=20
      NODE_MODULE_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/node_modules.zip
      BIO_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/bio.zip
      export JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/jpg.zip
      export JPG_NUMBER=508
      PCAP_DATA_FILES=1
  else
      export LOG_DATA_FILES=84
      export WAV_DATA_FILES=120
      NODE_MODULE_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/node_modules.zip
      BIO_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/bio.zip
      export JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/jpg.zip
      export JPG_NUMBER=1624
      PCAP_DATA_FILES=15
  fi
  
  ## WAV
  execute_step wav_step_1 wav_step_1_done_check "WAV zip download"
  execute_step wav_step_2 wav_step_2_done_check "WAV file generation"

  ## JPG
  execute_step jpg_step jpg_step_done_check "JPG Downloading"
  
  ## nginx logs
  execute_step nginx_logs_step_1 nginx_logs_step_1_done_check "NGINX logs Downloading"
  execute_step nginx_logs_step_2 nginx_logs_step_2_done_check "NGINX logs generated"
  
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
