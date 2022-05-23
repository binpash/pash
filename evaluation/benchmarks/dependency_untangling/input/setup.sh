#!/bin/bash

# exit when any command fails
#set -e

IN=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/
OUT=$PASH_TOP/evaluation/benchmarks/dependency_untangling/output/
IN_NAME=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/100G.txt

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

##
## This function checks if all the files in the arguments exist
## It returns 0 if all files exist, or 1 otherwise
##
files_exist_done_check()
{
    for file in "$@"; do
        if [ ! -f "$file" ]; then
            return 1
        fi
    done
    return 0
}

##
## This function executes a single idempotent step only if its check fails
##
## Requirements:
## - The step needs to be idempotent
## - The check needs to also check file sizes if there is concern of non-idempotence or failed download
##
execute_step()
{
    local step_fun=$1
    local step_done_check_fun=$2
    local step_desc=${3:-"Execution step"}

    # shellcheck disable=SC2086
    if ! eval $step_done_check_fun; then
        echo "$step_desc is not done, executing..."
        # shellcheck disable=SC2086
        eval $step_fun
        # shellcheck disable=SC2086
        eval $step_done_check_fun || { echo "ERROR: $step_desc failed!"; exit 1; }
    fi
    echo "$step_desc completed."
}

## Issues:
##
## - An overarching problem is that these take time in general, 
##   and therefore testing them out is not really feasible.
## - Another problem is that by doing that manually, 
##   we cannot get completely fine-grained. For example, we could
##   only copy the missing file _a la_ Rattle, instead of running
##   the whole step.
## - Another problem is that idempotence checking is hard to do manually.
## - Another issue is that generating the checks is cumbersome and error-prone.
## 

## Q: Can these checks be generated automatically?
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
    ## wget 
    wget -O wav.zip http://pac-n4.csail.mit.edu:81/pash_data/wav.zip
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


setup_dataset() {
  if [ "$1" == "--small" ]; then
      LOG_DATA_FILES=6
      export WAV_DATA_FILES=2
      NODE_MODULE_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/node_modules.zip
      BIO_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/bio.zip
      JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/small/jpg.zip
      PCAP_DATA_FILES=1
  else
      LOG_DATA_FILES=84
      export WAV_DATA_FILES=120
      NODE_MODULE_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/node_modules.zip
      BIO_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/bio.zip
      JPG_DATA_LINK=http://pac-n4.csail.mit.edu:81/pash_data/full/jpg.zip
      PCAP_DATA_FILES=15
  fi
  
  ## Step 1
  execute_step wav_step_1 wav_step_1_done_check "WAV zip download"

  ## Step 2
  execute_step wav_step_2 wav_step_2_done_check "WAV file generation download"

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
