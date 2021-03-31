INPUT=$PWD/input 
if [[ "$1" == "-c" ]]; then
    rm -rf input
    rm -rf output
    exit 
fi
# install dependencies
pkgs='ffmpeg unrtf imagemagick'
if ! dpkg -s $pkgs >/dev/null 2>&1; then
  sudo apt-get install $pkgs
fi

# create directories
mkdir -p output
mkdir -p $INPUT
if [ ! -f $INPUT/tomhannen-20080409.tgz ]; then
    cd $INPUT
    echo "Fetching Dataset"
    wget http://www.repository.voxforge1.org/downloads/SpeechCorpus/Trunk/Audio/Original/48kHz_16bit/tomhannen-20080409.tgz
    tar xf tomhannen-20080409.tgz
    #rm -rf wav mp3
    mkdir wav rtf
    cd tomhannen-20080409/wav
    # total 5.7 size of audio files
    for i in *.wav; do
        FILE=$(basename "$i")
        for x in {1..20}; do cp $i "../../wav/$i$x.wav"; done
    done
    #cleanup
    rm -rf ../../tomhannen-20080409
    
fi
# fetch sample.rtf
if [ ! -f $INPUT/sample.rtf ]; then
    cd $INPUT
    wget https://jeroen.github.io/files/sample.rtf 
    for i in {0..10000}
    do
        cp sample.rtf rtf/sample$i.rtf
    done
    
fi

# generate images
if [ ! -f $INPUT/jpg1.tar.gz ]; then
    cd $INPUT
    echo "Fetching Dataset"                                                   
    #wget ftp://ftp.inrialpes.fr/pub/lear/douze/data/jpg1.tar.gz
    tar xf jpg1.tar.gz
    rm -rf tmp1 
    mkdir -p tmp1
    cd jpg/
    for filename in *.jpg; do                                                     
        cp $filename ../tmp1/${filename}_copy.jpg                                  
    done
    cd ../tmp1
    for filename in *.jpg; do                                                     
        cp $filename ../jpg/
    done                                                                          
    echo "JPGs copied"
    rm -rf ../tmp1
fi
# gen dummy files
if [ ! -d $INPUT/tmp ]; then
    cd $INPUT
    mkdir tmp
    seq -w 1 100000 | xargs -P 100 -I{} sh -c 'num=$(echo {} | sed 's/^0*//');val=$(($num % 1000)); touch tmp/f{}.$val;'
    echo "Generated 100000 empty files"
fi

me=`basename "$0"`
echo "$me completed"
