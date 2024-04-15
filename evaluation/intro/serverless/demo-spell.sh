#!/bin/sh

cd "$(dirname $0)"

S3_FILE="input/200M.txt"
FILE="/tmp/200M.txt"
S3_DICT="input/sorted_words"
DICT="/tmp/sorted_words"

python3 aws/s3-get-object.py $S3_FILE $FILE
python3 aws/s3-get-object.py $S3_DICT $DICT

cat "$FILE" | tr A-Z a-z | tr -cs A-Za-z '\n' | sort | uniq | comm -13 $DICT - >/dev/null
