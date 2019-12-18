#!/bin/bash

# https://www.biostars.org/p/43677/
# https://github.com/h3abionet/h3agatk
# https://docs.google.com/document/d/1siCZrequI4plggz3ho351NnX57CoyCJl9GWp3azlxfU/edit#
bwa mem -M -p -t [num_threads] \
      -R "@RG\tID:1\tPL:ILLUMINA\tPU:pu\tLB:group1\tSM:SAMPLEID" \
      [reference_fasta] \
      [input_fastq] > [output]

bwa mem genome.fa reads.fastq | samtools sort -o output.bam -

# https://www.biostars.org/p/43677/
bwa aln -t 4 ./hg19.fasta ./s1_1.fastq > ./s1_1.sai
bwa aln -t 4 ./hg19.fasta ./s1_2.fastq > ./s1_2.sai
bwa sampe ./hg19.fasta ./s1_1.sai ./s1_2.sai ./s1_1.fastq ./s1_2.fastq |
  samtools view -Shu - |
  samtools sort - - |
  samtools rmdup -s - - |
  tee s1_sorted_nodup.bam |
  bamToBed > s1_sorted_nodup.bed

# 4 cores, -M is for Picard compatibility
bwa mem -M -t 4 ./hg19.fasta ./s1_1.fastq ./s1_2.fastq > s1.sam

samtools merge - *.bam |
# tee merged.bam |
  samtools rmdup - - |
# tee rmdup.bam |
  samtools mpileup - uf ./hg19.fasta - |
  bcftools view -bvcg - | gzip > var.raw.bcf.gz

bwa sampe ./hg19.fasta <(bwa aln -t 4 ./hg19.fasta ./s1_1.fastq) <(bwa aln -t 4 ./hg19.fasta ./s1_2.fastq) ./s1_1.fastq ./s1_2.fastq | samtools view -Shb /dev/stdin > s1.bam
