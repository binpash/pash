# convert all rtf to txt
rtf_to_txt()
(
    find $1 -name "*.rtf" | xargs -I {} unrtf {} --text > /dev/null
)


rtf_to_txt  input/rtf
