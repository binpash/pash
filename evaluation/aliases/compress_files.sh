# compress all the files in a directory using dd and tar
INPUT=${INPUT:-$PASH_TOP/evaluation/aliases/input/rtf}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/aliases/output}
cd $INPUT
# get all rtf and compress them
find . -name "*.rtf" | xargs -P16 -I {} sh -c "dd if={} bs=1 status=none > '{}f'; tar -zcf {}.tar.gz {}f; rm {}f; mv {}.tar.gz $OUTPUT" sh {}
