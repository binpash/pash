DIR=$PASH_TOP/evaluation/scripts/input/
mkdir -p $DIR
cd $DIR
if [[ ! -f nginx.zip ]]
then
    wget -O nginx.zip https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/3QBYB5/NXKB6J
    unzip nginx.zip 
fi

