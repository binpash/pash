INPUT=${INPUT:-$PASH_TOP/evaluation/bio/input/}
# remove adapter
find ${INPUT} -name "*.fastq" |  sort | uniq | xargs -I {} cutadapt -a AGATCGGAAGAGCACAC {} >  /dev/null
