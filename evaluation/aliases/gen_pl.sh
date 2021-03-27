# generate a playlist
INPUT=${INPUT:-$PASH_TOP/evaluation/aliases/input/}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/aliases/output}
find $1 -type f -name *.mp3 -o -name *.wav | sort > $OUTPUT/playlist.pls
