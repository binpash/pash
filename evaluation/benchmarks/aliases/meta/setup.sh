#!/bin/bash

# exit when any command fails
set -e


IN=$PASH_TOP/evaluation/benchmarks/aliases/meta/
OUT=$PASH_TOP/evaluation/benchmarks/aliases/meta/out
if [[ "$1" == "-c" ]]; then
    rm -rf ${IN}/jpg
    rm -rf ${IN}/rtf
    rm -rf ${IN}/wav
    rm -rf ${OUT}
    exit 
fi




mkdir -p ${OUT} 
cd $IN

# mp3 dataset
# many, small files
# if [ ! -f tomhannen-20080409.tgz ]; then                                                 
#   echo "Fetching Dataset"                                                   
#   wget http://www.repository.voxforge1.org/downloads/SpeechCorpus/Trunk/Audio/Original/48kHz_16bit/tomhannen-20080409.tgz
#   tar xf tomhannen-20080409.tgz
# fi
# total 5.7 size of audio files

mkdir -p wav
cd wav

if [[ ! -f file_example_WAV_1MG.wav ]]; then
  wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_1MG.wav
  wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_2MG.wav
  wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_5MG.wav
  wget https://file-examples-com.github.io/uploads/2017/11/file_example_WAV_10MG.wav

  for f in *.wav; do
    FILE=$(basename "$f")
    for i in {1..20}; do 
      echo copying to $f$i.wav
      cp $f $f$i.wav
    done
  done
fi

cd  $IN
mkdir -p rtf
cd rtf
if [[ ! -f sample.rtf ]]; then
  wget https://jeroen.github.io/files/sample.rtf
  for i in {0..10000}; do
    cp sample.rtf $i.rtf
  done
fi

cd  $IN
mkdir -p jpg
cd jpg
if [ ! -f 107801.jpg ]; then
  curl ftp://ftp.inrialpes.fr/pub/lear/douze/data/jpg1.tar.gz | tar -xz
  mv jpg/* . && rm -r jpg/
fi                                                                            

exit 

# mkdir -p tmp                                                                  
# cd jpg/
# for filename in *.jpg; do                                                     
#     cp $filename ../tmp/${filename}_copy.jpg                                  
# done
# 
# cd ../tmp
# for filename in *.jpg; do                                                     
#     cp $filename ../jpg/
# done                                                                          
# echo "JPGs copied"
# rm -rf ../tmp
# cd ..
# rm -rf tmp
# 
# mkdir tmp
# seq -w 1 100000 | xargs -P 100 -I{} sh -c 'num=$(echo {} | sed 's/^0*//');val=$(($num % 1000)); touch tmp/f{}.$val;'
# echo "Generated 100000 empty files"
