# Here are sample steps to generate a single paired read from hg19:
# https://www.biostars.org/p/150010/
INPUT=${INPUT:-$PASH_TOP/evaluation/bio/input/}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/bio/output}
cd ${INPUT}
# filter out a single chromosome and index it, e.g.
samtools faidx ${INPUT}/human_g1k_v37.fasta 20 > ${OUTPUT}/human_g1k_v37_chr20.fasta
bowtie2-build ${OUTPUT}/human_g1k_v37_chr20.fasta ${OUTPUT}/homo_chr20
#simulate a single read sample, e.g. here is for a single (-N 1) paired read:
${INPUT}/wgsim/wgsim -N 1 ${OUTPUT}/human_g1k_v37_chr20.fasta ${OUTPUT}/single.read1.fq ${OUTPUT}/single.read2.fq > ${OUTPUT}/wgsim.out
#generate the sam, e.g.
bowtie2 -x ${OUTPUT}/homo_chr20 -1 ${OUTPUT}/single.read1.fq -2 ${OUTPUT}/single.read2.fq -S ${OUTPUT}/single_pair.sam
#generate a bam
samtools view -b -S -o ${OUTPUT}/single_pair.bam ${OUTPUT}/single_pair.sam 
#sort and index it
samtools sort ${OUTPUT}/single_pair.bam -o ${OUTPUT}/single_pair.sorted.bam
# this seems to not affect the file, but in other cases, its indeed needed
samtools index ${OUTPUT}/single_pair.sorted.bam 

