INPUT=${INPUT:-$PASH_TOP/evaluation/bio/input/bam}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/bio/output}
cd ${INPUT}
find . -name "*.bam" | xargs -I {} samtools view -h -o ${OUTPUT} {}
