# convert wat to mp3
INPUT=${INPUT:-$PASH_TOP/evaluation/aliases/input/wav}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/aliases/output}
cd $INPUT
find . -name "*.wav" | xargs -I {} -P 16 ffmpeg -y -hide_banner -loglevel error -i {} -ab 192000 $OUTPUT/{}.mp3
