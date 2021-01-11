mkfifo {1,2,3}grams

bigrams_aux()
{
    ( mkfifo s2 > /dev/null ) ;
    ( mkfifo s3 > /dev/null ) ;

    sed '$d' s2 > s3 &
    tee s2 |
        tail +2 |
        paste s3 -
    rm s2
    rm s3
}

bigram_aux_map()
{
    IN=$1
    OUT=$2
    AUX_HEAD=$3
    AUX_TAIL=$4

    s2=$(mktemp -u)
    aux1=$(mktemp -u)
    aux2=$(mktemp -u)
    aux3=$(mktemp -u)
    temp=$(mktemp -u)

    mkfifo $s2
    mkfifo $aux1
    mkfifo $aux2
    mkfifo $aux3

    ## New way of doing it using an intermediate file. This is slow
    ## but doesn't deadlock
    cat $IN > $temp

    sed '$d' $temp > $aux3 &
    cat $temp | head -n 1 > $AUX_HEAD &
    cat $temp | tail -n 1 > $AUX_TAIL &
    cat $temp | tail +2 | paste $aux3 - > $OUT &

    wait

    rm $temp
    rm $s2
    rm $aux1
    rm $aux2
    rm $aux3
}

bigram_aux_reduce()
{
    IN1=$1
    AUX_HEAD1=$2
    AUX_TAIL1=$3
    IN2=$4
    AUX_HEAD2=$5
    AUX_TAIL2=$6
    OUT=$7
    AUX_HEAD_OUT=$8
    AUX_TAIL_OUT=$9

    temp=$(mktemp -u)

    mkfifo $temp

    cat $AUX_HEAD1 > $AUX_HEAD_OUT &
    cat $AUX_TAIL2 > $AUX_TAIL_OUT &
    paste $AUX_TAIL1 $AUX_HEAD2 > $temp &
    cat $IN1 $temp $IN2 > $OUT &

    wait

    rm $temp
}


trigrams_aux()
{
    s2=$(mktemp -u)
    s3=$(mktemp -u)

    mkfifo $s2 $s3

    tee $s2 |
        tail +2 |
        paste $s2 - |
        tee $s3 |
        cut -f 1 |
        tail +3 |
        paste $s3 - |
        sed "\$d" |
        sed "\$d"

    rm $s2 $s3
}


extract_text()
{
    while read -r line
    do
        cat $line |
            iconv -c -t ascii//TRANSLIT |
            pandoc +RTS -K64m -RTS --from html --to plain --quiet
    done
}


cat $IN |
  sed "s#^#$WIKI#" |
  extract_text |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf $WEB_INDEX_DIR/stopwords.txt |
  $WEB_INDEX_DIR/stem-words.js |
  tee 3grams 2grams 1grams > /dev/null &

cat 1grams |
    sort |
    uniq -c |
    sort -rn > 1-grams.txt &

cat 2grams |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    bigrams_aux |
    sort |
    uniq -c |
    sort -rn > 2-grams.txt &

cat 3grams |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    trigrams_aux |
    sort |
    uniq -c |
    sort -rn # > 3-grams.txt

rm {1,2,3}grams
