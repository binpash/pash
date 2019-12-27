#!/bin/bash

# Find differences between two genome sequences---a a paired Illumina sequencing
# read  (FASTQ files)  and an  assembled  reference genome  from GenBank  (e.g.,
# Pasteurella multocida). The reads are aligned to the reference,  and sorted by
# coordinate. Instead of saving the BAM file, we pipe it directly to a series of
# BCF tool  steps. Note the use  of -l 0  and -Ou to  keep the piped data  in an
# uncompressed  form, to  avoid  repeated  compression/decompression steps.  The
# --min-MQ 60 ensures only uniquely mapped reads are used. The final filter step
# removes  low  quality  variant  calls,  heterozygous  calls  (this  is  haploid
# bacteria), and any regions with less than 10 supporting reads.

# Requires: samtools, minimap2, bcftools
# Data: http://ndr.md/data/bio/R1.fastq.gz http://ndr.md/data/bio/R2.fastq.gz  http://ndr.md/data/bio/ref.fa.gz

# https://github.com/samtools/samtools/releases/latest
# https://github.com/lh3/minimap2
# http://thegenomefactory.blogspot.com/2018/10/a-unix-one-liner-to-call-bacterial.html

CPUS=4
REF=ref.fa
R1=R1.fastq.gz
R2=R2.fastq.gz

minimap2 -a -x sr -t "$CPUS" "$REF" "$R1" "$R2" |             # align reads to the reference
  samtools sort -l 0 --threads "$CPUS" |                      # sort reads by coordinate
  bcftools mpileup -Ou -B --min-MQ 60 -f "$REF" - |           # multi-way pileup producing genotype likelihoods
  bcftools call -Ou -v -m - |                                 # SNP/indel calling 
  bcftools norm -Ou -f "$REF" -d all - |                      # left-align and normalize indels
  bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' |    # removes  low-quality  variant  calls, etc
  bcftools stats |                                            # produce VCF/BCF stats
  grep '^SN' |                                                # 
  cut -f3-
