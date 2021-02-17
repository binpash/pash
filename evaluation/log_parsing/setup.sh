DIR=$PASH_TOP/evaluation/scripts/input/log_parsing
mkdir -p $DIR
cd $DIR
if [[ ! -f 201011271400.dump.gz ]]
then
    wget http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201011271400.dump.gz
fi
gunzip 201011271400.dump.gz
