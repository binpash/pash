#!/bin/bash

# http://thegenomefactory.blogspot.com/2018/10/a-unix-one-liner-to-call-bacterial.html
CPUS=4
REF=ref.fa
R1=R1.fastq.gz
R2=R2.fastq.gz

minimap2 -a -x sr -t "$CPUS" "$REF" "$R1" "$R2" \
 | samtools sort -l 0 --threads "$CPUS" \
 | bcftools mpileup -Ou -B --min-MQ 60 -f "$REF" - \
 | bcftools call -Ou -v -m - \
 | bcftools norm -Ou -f "$REF" -d all - \
 | bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' 
 > variants.vcf

 # Then:
# bcftools stats variants.vcf | grep '^SN' | cut -f3-
