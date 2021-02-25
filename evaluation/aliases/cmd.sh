convert_to_mp3()
(
   cd $1
   find . -name "*.wav" | xargs -I {} -P 16 ffmpeg -i {} -ab 192000 ../mp3/{}.mp3
   cd ../mp3
   # sort based on name
    for file in *.mp3; do
      dir=${file:0:2}
      mkdir -p "$dir"  &&      mv -iv "$file" "./$dir"
    done


    cd - 
)
#convert_to_mp3 $1
# generate a playlist
gen_playlist()
(
  cd $1
  find . -type f -name *.mp3 -o -name *.wav | sort > playlist.pls
)

grep_words_ing()
(
    # grep words ending in ing
    find . -name "*.txt" | xargs -n 1 cat | tr -sc 'A-Za-z' '\n' | tr 'A-Z' 'a-z' | grep 'ing$' | sort | uniq -c | sort -n -r 
)

delete_movies()
(   
    # retrieve files with an extension, move them to a seperate folder
    mkdir -p ~/movies_to_delete && find  -type f -regex ".*\.\(mkv\|avi\|mov\|mpg\|wav\)"  -print0 | xargs -0 -I{} dirname {}  |xargs -I{} mv {} ~/movies_to_delete
)

convert_to_mp3 "$PASH_TOP/evaluation/scripts/input/aliases/wav"
gen_playlist "$PASH_TOP/evaluation/scripts/input/aliases/"
