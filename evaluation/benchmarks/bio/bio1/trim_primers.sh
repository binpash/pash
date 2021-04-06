# trim primers
INPUT=${INPUT:-$PASH_TOP/evaluation/bio/input/}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/bio/output}
cd ${INPUT}
find . -maxdepth 1 -name "*.fastq" | xargs -I {}  cutadapt -a TCCTCCGCTTATTGATAGC -o ${OUTPUT}/{}\_trimmed.fastq {}; 

