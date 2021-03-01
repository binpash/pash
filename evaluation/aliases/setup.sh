# remove this line since it hangs timeout 2 -i brightness | cut -f2 -d "" "" > ~/.currentbrightness; cat\
# execute once
#sed -e '46106d' $1                                                               
# likely-longest-pipelines.txt is locally stored
# strip the number of pipe column
cut -f1 -d " " --complement  likely-longest-pipelines.txt > $PASH_TOP/evaluation/scripts/input/generated.file

## mp3 dataset ##
apt-get install ffmpeg unrtf imagemagick
PW=$PASH_TOP/evaluation/scripts/input/aliases
mkdir -p $PW
cd $PW
wget http://www.repository.voxforge1.org/downloads/SpeechCorpus/Trunk/Audio/Original/48kHz_16bit/tomhannen-20080409.tgz
tar xf tomhannen-20080409.tgz
rm -rf wav mp3
mkdir wav mp3 rtf
cd tomhannen-20080409/wav
# total 5.7 size of audio files
for i in *.wav; do
  FILE=$(basename "$i")
  for x in {1..2}; do cp $i "../../wav/$i$x.wav"; done
done
rm -rf ../../tomhannen-20080409
cd  $PW
# fetch all gutenberg books here ?
# NLP grep
rm *.tar.gz
wget -O nlp.tar.gz https://wordnetcode.princeton.edu/3.0/WNdb-3.0.tar.gz
tar xf nlp.tar.gz

wget https://jeroen.github.io/files/sample.rtf
for i in {0..10000}
do
    cp sample.rtf rtf/sample$i.rtf
done

if [ ! -f jpg1.tar.gz ]; then                                                 
    echo "Fetching Dataset"                                                   
    wget -O jpg1.tar.gz ftp://ftp.inrialpes.fr/pub/lear/douze/data/jpg1.tar.gz
fi                                                                            
tar xf jpg1.tar.gz                                                            
mkdir -p tmp                                                                  
cd jpg/                                                                       
for filename in *.jpg; do                                                     
    echo $filename                                                            
    cp $filename ../tmp/${filename}_copy.jpg                                  
done                                                                          

mv ../tmp/* .
rm -rf ../tmp
cd ..
