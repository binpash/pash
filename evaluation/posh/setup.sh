mkdir -p $PASH_TOP/evaluation/scripts/input/posh
cd $PASH_TOP/evaluation/scripts/input/posh/
rm -rf output
rm -rf jpg
mkdir -p output
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
wget -O log1.zip http://www.sec.gov/dera/data/Public-EDGAR-log-file-data/2017/Qtr1/log20170314.zip 
wget -O log2.zip http://www.sec.gov/dera/data/Public-EDGAR-log-file-data/2017/Qtr1/log20170315.zip
unzip log1.zip
rm README.txt
unzip log2.zip

