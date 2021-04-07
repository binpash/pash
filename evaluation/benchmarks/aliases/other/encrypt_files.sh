# compress and encrypt all files in a directory 
INPUT=${INPUT:-$PASH_TOP/evaluation/aliases/input/rtf}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/aliases/output}
cd $INPUT
find . -name "*.rtf" | xargs -I {} sh -c "tar -czf - {} | openssl enc -e -pbkdf2 -out {}.enc; mv {}.enc $OUTPUT" sh {}
