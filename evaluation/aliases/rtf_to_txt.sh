INPUT=${INPUT:-$PASH_TOP/evaluation/aliases/input/rtf}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/aliases/output}
# convert all rtf to txt
find $INPUT -name "*.rtf" | xargs -I {} unrtf {} --text > /dev/null
