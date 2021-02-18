DIR=$PASH_TOP/evaluation/scripts/input/
mkdir -p $DIR
cd $DIR
if [[ ! -f 201011271400.dump.gz ]]
then
    # too slow
    wget http://mawi.wide.ad.jp/mawi/samplepoint-F/2010/201011271400.dump.gz
    gunzip 201011271400.dump.gz
fi

