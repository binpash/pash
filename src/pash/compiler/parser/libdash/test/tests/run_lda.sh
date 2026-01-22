#!/bin/bash

export PYTHONIOENCODING=utf8

if test $# -ne 0;
then
    KS="$*";
else
    KS="50 75 100 125 150 175 200"
fi

DIR=`date "+%Y-%m-%d_%H:%M"`
START=`date "+%Y-%m-%d %H:%M"`

# TODO error handling

echo "SETTING UP"
mkdir ${DIR}

echo "PARSING"
python parse.py

for dat in abstracts.dat vocab.dat docs.dat; do
    mv ${dat} ${DIR}
done

# we don't want to lose this one!
cp stopwords.dat ${DIR}

echo "RUNNING LDA"

ABS=${DIR}/abstracts.dat

for k in ${KS}; do
    lda est 1/50 ${k} settings.txt ${ABS} seeded ${DIR}/lda${k} &
    echo lda${k} >>${DIR}/.gitignore
done

wait
echo "PROCESSING TOPICS"

for k in ${KS}; do
    python debug_topics.py ${DIR} ${k} > ${DIR}/lda${k}_topics.txt
done

echo "GENERATING CSV"

for i in ${DIR}/lda*; do
    test -d ${i} && python post.py ${i}/final.gamma ${DIR}/docs.dat > ${i}.csv
    test -d ${i} && python by_year.py ${i}/final.gamma ${DIR}/docs.dat > ${i}_by_year.csv
done

echo "MOVING TO OUTPUT DIRECTORY"
mv ${DIR} ../out

echo "DONE"
echo All done. Started at ${START}, done at `date "+%Y-%m-%d %H:%M"`.
