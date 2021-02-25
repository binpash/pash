# remove this line since it hangs timeout 2 -i brightness | cut -f2 -d "" "" > ~/.currentbrightness; cat\
# execute once
#sed -e '46106d' $1                                                               
# likely-longest-pipelines.txt is locally stored
# strip the number of pipe column
cut -f1 -d" " --complement  likely-longest-pipelines.txt > $PASH_TOP/evaluation/scripts/input/generated.file

## mp3 dataset ##
apt-get install ffmpeg
PW=$PASH_TOP/evaluation/scripts/input/aliases
mkdir -p $PW
cd $PW
wget http://www.repository.voxforge1.org/downloads/SpeechCorpus/Trunk/Audio/Original/48kHz_16bit/tomhannen-20080409.tgz
tar xf tomhannen-20080409.tgz
rm -rf wav mp3
mkdir wav mp3
cd tomhannen-20080409/wav
# total 5.7 size of audio files
for i in *.wav; do
  FILE=$(basename "$i")
  for x in {1..200}; do cp $i "../../wav/$i$x.wav"; done
done
rm -rf ../../tomhannen-20080409
