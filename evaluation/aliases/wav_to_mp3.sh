convert_to_mp3()
(
    cd $1
   find . -name "*.wav" | xargs -I {} -P 16 ffmpeg -y -hide_banner -loglevel error -i {} -ab 192000 $2/{}.mp3
)

convert_to_mp3 input/wav $PWD/output
