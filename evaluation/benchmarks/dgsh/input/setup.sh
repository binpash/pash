if [[ $1 == "-c" ]]; then
    rm -rf input/
    exit 1
fi
mkdir -p input
cd input
if [[ ! -f dblp.xml ]]; then
    wget https://dblp.uni-trier.de/xml/dblp.xml.gz
    gunzip dblp.xml.gz
    cat dblp.xml | head -n 35 > mini.xml
fi

if [[ ! -f fid ]]; then
    wget http://www.bmrb.wisc.edu/ftp/pub/bmrb/timedomain/bmr6443/timedomain_data/c13-hsqc/june11-se-6426-CA.fid/fid
    wget https://www.stats.govt.nz/assets/Uploads/International-trade/International-trade-December-2020-quarter/Download-data/international-trade-december-2020-quarter-csv.zip
    unzip *.zip
    mv final/* .
    rm -rf final
fi



