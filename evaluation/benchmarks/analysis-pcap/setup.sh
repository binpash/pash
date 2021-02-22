DIR=$PASH_TOP/evaluation/scripts/input/
mkdir -p $DIR
cd $DIR
if [[ ! -f 201011271400.dump.gz ]]
then
    # too slow
    wget http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201011271400.dump.gz
    gunzip 201011271400.dump.gz
    wget https://mcfp.felk.cvut.cz/publicDatasets/IoT-23-Dataset/IndividualScenarios/CTU-IoT-Malware-Capture-7-1/2018-07-20-17-31-20-192.168.100.108.pcap
fi

