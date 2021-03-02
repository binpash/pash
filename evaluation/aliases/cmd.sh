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

# we can maybe move this to another script ?
#grep_words()
#(
#    # the first n books
#    dict_books=1000
#    # our dictionary file
#    my_dict=$1/my_dict
#    find $1 -name "*.txt" | head -n $dict_books > $1/${dict_books}.entries
#    ################################################
#    # Generate dictionary based on all the files 
#    ################################################
#    cat $1/${dict_books}.entries | xargs -n 1 cat | tr -sc 'A-Za-z' '\n' | tr 'A-Z' \
#    'a-z' | sort | uniq > $my_dict
#    echo "Created dictionary from allthe books"
#    #################################################
#    # grep words ending in ing
#    #################################################
#    cat $my_dict | grep 'ing$' | sort | uniq -c | sort -n -r  > $1/res
#    echo "Found all the -ings"
#    #################################################
#    # find all duplicate words across multiple files
#    #################################################
#    file_expr="$1/*.txt"; echo $file_expr > $1/f; sort $1/f | sed 's/^\s*//; s/\s*$//; /^\s*$/d' | \
#  uniq -d | while read dup_line; do grep -Hn "^\s*$dup_line\s*$" $file_expr; \
#done| sort -t: -k3 -k1,2 | awk -F: '{ file=$1; line=$2; $1=$2="";gsub(/(^[\t]+)| \
#([ \t]+$)/,"",$0); if (prev != "" && prev != $0) printf ("\n"); printf("\033[0;33m%s  \
#    (line %s)\033[0m: %s\n", file, line, $0); prev=$0; }' > $1/out
#    echo "Found all the duplicates"
#    #################################################
#    # Get words that don't exist in the dictionary
#    ################################################
#    # get all the words from the files that don't exist in dict 
#    #################################################
#    comm -1 $my_dict /usr/share/dict/words > $1/words
#    echo "Found words uniq to our dict"
#    ##################################################
#    # NLP ?
#    ##################################################
#    adj=$1/conv.data.adj
#    adv=$1/conv.data.adv
#    noun=$1/conv.data.noun
#    verb=$1/conv.data.verb
#    egrep -o "^[0-9]{8}\s[0-9]{2}\s[a-z]\s[0-9]{2}\s[a-zA-Z]*\s" $1/dict/data.adj | cut \
#    -d ' ' -f 5 > $adj
#    egrep -o "^[0-9]{8}\s[0-9]{2}\s[a-z]\s[0-9]{2}\s[a-zA-Z]*\s" $1/dict/data.adv | cut \
#    -d ' ' -f 5 > $adv
#    egrep -o "^[0-9]{8}\s[0-9]{2}\s[a-z]\s[0-9]{2}\s[a-zA-Z]*\s" $1/dict/data.noun | cut \
#    -d ' ' -f 5 > $noun
#    egrep -o "^[0-9]{8}\s[0-9]{2}\s[a-z]\s[0-9]{2}\s[a-zA-Z]*\s" $1/dict/data.verb | cut \
#    -d ' ' -f 5 > $verb
#    comm -12 $(cat $adj  | sort > $1/.tmp; echo $1/.tmp) $my_dict >$1/r1
#    comm -12 $(cat $adv  | sort > $1/.tmp; echo $1/.tmp) $my_dict >$1/r2
#    comm -12 $(cat $noun | sort > $1/.tmp; echo $1/.tmp) $my_dict >$1/r3
#    comm -12 $(cat $verb | sort > $1/.tmp; echo $1/.tmp) $my_dict >$1/r4
#    # Get verbs from the global dict
#    # get the verbs from each file that exists in global dict
#    # count the number of uniq verbs for each file
#    cat $1/${dict_books}.entries | xargs -I {} sh -c 'comm -12 $(cat "$1" | \
#      tr -sc "A-Za-z" "\n" | tr "A-Z" "a-z" | sort | uniq > "$1".lol; \
#      echo "$1".lol) $PASH_TOP/evaluation/scripts/input/aliases/r4 | wc -l > "$1".wc' sh {} 
#    echo "Generated word count files"
#    
#
#)

delete_movies()
(   
    # retrieve files with an extension, move them to a seperate folder
    mkdir -p ~/movies_to_delete && find  -type f -regex ".*\.\(mkv\|avi\|mov\|mpg\|wav\)"  -print0 | xargs -0 -I{} dirname {}  |xargs -I{} mv {} ~/movies_to_delete
)

# convert to jpg to psd
# can we optimize this ?
jpg_to_psd()
(
    cd $1
    mkdir converted;
    find . -name "*.jpg" | xargs -P 16 -I {} sh -c "convert {} converted/$(echo {} | sed -e 's/\.jpg/\.psd/g') &&  echo Converted {} to $(echo {} | sed -e 's/\.jpg/\.psd/g') ...;" sh {}
)

# Resize image 
resize_image()
(
    find $1 -name *.jpg -o -name *.png -o -name *.jpeg -o -name *.JPG -o -name *.PNG \
        -o -name *.JPEG | xargs -P 16 -I {} convert -resize 70% {} {}.70
)

# convert all rtf to txt
rtf_to_txt()
(
    cd $1
    find . -name "*.rtf" | xargs -I {} unrtf {} --text > /dev/null
)

# count how many times each file type exist in a directory
# 1GB files take about 14s on i7 8th gen
get_type_count()
(
    cd $1
    find . -type f | while read f; do echo ""${f##*.}""; done | sed ""/^\s*$/d"" | sort | uniq -c | sort -rn > /dev/null
)

# calculate a hash? can we change it to calculate hashes for all the files?
get_hash()
(
    head -c32 /dev/urandom | openssl dgst -sha256 -binary -hmac $(xxd -p -l32 -c32 /dev/urandom) | base64 | cut -b-32
)

# compress all the files in a directory using dd and tar
compress_files()
(
    cd $1
    find . -name "*.rtf" | xargs -I {} sh -c "dd if={} bs=1 skip=24 > '{}f'; tar -zcvf {}.tar.gz {}f; rm {}f" sh {}
)

# compress and encrypt all files in a directory 
encrypt_files()
(
    cd $1
    find . -name "*.rtf" | xargs -I {} sh -c "tar -czf - {} | openssl enc -e -pbkdf2 -out {}.enc" sh {}
)

# no idea what this does, but it takes time
func()
(
while read line; 
do 
    count=1; 
    out=$(echo "$line"|tr "|" "\n");  
    echo $out |
        while IFS= read -r i; 
        do echo "$count $i"; 
            count=$((count+1)); 
        done; 
        echo "-----------------------------------------------      -----------------------------------------------------"; 
    done
)

#convert_to_mp3 "$PASH_TOP/evaluation/scripts/input/aliases/wav"
#gen_playlist "$PASH_TOP/evaluation/scripts/input/aliases/"
#resize_image "$PASH_TOP/evaluation/scripts/input/aliases/jpg"
#jpg_to_psd "$PASH_TOP/evaluation/scripts/input/aliases/jpg"
#rtf_to_txt  "$PASH_TOP/evaluation/scripts/input/aliases/rtf"
#get_type_count "$PASH_TOP/evaluation/scripts/input/aliases/tmp"
#compress_files $PASH_TOP/evaluation/scripts/input/aliases/rtf
#encrypt_files $PASH_TOP/evaluation/scripts/input/aliases/rtf

