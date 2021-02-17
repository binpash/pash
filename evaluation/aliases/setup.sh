# remove this line since it hangs timeout 2 -i brightness | cut -f2 -d "" "" > ~/.currentbrightness; cat\
# execute once
#sed -e '46106d' $1                                                               
#rm -f $PASH_TOP/evaluation/scripts/input/aliases/*
mkdir -p $PASH_TOP/evaluation/scripts/input/aliases/
# likely-longest-pipelines.txt is locally stored
# strip the first column
cut -f1 -d" " --complement  likely-longest-pipelines.txt > $PASH_TOP/evaluation/scripts/input/aliases/generated.file
