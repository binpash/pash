# convert fastq to fasta format
# It recognizes the extension .fasta and it converts the input to fasta.gz format
INPUT=${INPUT:-$PASH_TOP/evaluation/bio/input/}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/bio/output}
cd ${INPUT}
find . -maxdepth 1 -name "*.fastq" | xargs -I {} cutadapt -o ${OUTPUT}/{}.fasta.gz {}
