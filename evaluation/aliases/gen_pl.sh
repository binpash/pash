# convert_to_mp3 $1
# generate a playlist
gen_playlist()
(
    find $1 -type f -name *.mp3 -o -name *.wav | sort > $2
)

gen_playlist input/ $PWD/output/playlist.pls
