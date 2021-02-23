# dataset http://www.repository.voxforge1.org/downloads/SpeechCorpus/Trunk/Audio/Original/48kHz_16bit/tomhannen-20080409.tgz
convert_to_mp3()
(
    cd $1
    for FILE in *.wav ;
    do
        echo $FILE
        ffmpeg -i $FILE -f mp3 -ab 192000 ${FILE}.mp3 || break;
    done
)
#convert_to_mp3 $1
# generate a playlist
gen_playlist()
(
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
