#### Ported ####
# https://dfzljdn9uc3pi.cloudfront.net/2013/203/1/Supplement_S2.pdf
set -e
cd $PASH_TOP/evaluation/benchmarks/bio/bio1/input/
ls *.R1.fq > namelist 
sed -i 's/.R1.fq//g' namelist 
NAMES=( `cat "namelist" `)
mkdir -p assembly 
# Trims raw files two different ways.
# First way removes any reads with substantial amounts of adapter, but does no
# quality trimming. These reads are used for assembly and must be uniform lengths
# Second way removes adapters and does quality trimming. These reads will be
# used for mapping. 
for i in "${NAMES[@]}"
do
  echo $i
  Trim/trim_galore --paired -q 0 --length 90 -a   GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2   GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 20   ${i}.R1.fq ${i}.R2.fq --output_dir ./assembly
  Trim/trim_galore --paired -q 20 --length 20 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATATCGTATGCCGTCTTCTGCTTG -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG --stringency 10 $i.R1.fq $i.R2.fq
done 

# Renaming trimmed files to simpler names
for i in "${NAMES[@]}"
do
  mv $i.R1_val_1.fq $i.1.fq
  mv $i.R2_val_2.fq $i.2.fq
done

### Assembly ###
# These parameters could be further optimized for particular taxa
# First step concatenates reads into one forward and one reverse fastq file
cat ./assembly/*.R1_val_1.fq > forward
cat ./assembly/*.R2_val_2.fq > reverse
# Rainbow now clusters and assembles
rainbow/rainbow cluster -1 forward -2 reverse > cat.rbcluster.out 2> log
# we can add -f $1 but im not good with maths
rainbow/rainbow div -i cat.rbcluster.out -o cat.rbdiv.out 
rainbow/rainbow merge -a -i cat.rbdiv.out -o cat.rbasm.out -N 1000
perl rainbow/select_best_rbcontig.pl cat.rbasm.out > rainbowf
# Renames contigs to sequential numbers for simplicity
fastx_renamer -n COUNT -i rainbowf -o reference 
## Mapping
# Use BWA to index reference
bwa-0.7.17/bwa index -a bwtsw reference
# Use BWA to map reads to reference.
### These parameters could be further optimized for particular taxa
for i in "${NAMES[@]}"
do
  bwa-0.7.17/bwa mem reference $i.1.fq $i.2.fq -t 32 -a -T 10 > $i.sam
done 
#Convert Sam to Bam and remove low quality, ambiguous mapping
for i in "${NAMES[@]}"
do
  samtools view -bS -q15 $i.sam > $i.bam
  samtools sort $i.bam -o $i
done
# Index reference for SAMtools
samtools faidx reference
# sort the Sample1.bam cause it sucks. The file needs to be sorted in that
# way before index is called
samtools sort -m 2G -@ 4 Sample1.bam -o lala
mv lala Sample1.bam
# index the bamfile
samtools index Sample1.bam
samtools mpileup -D -f reference *.bam >mpileup 
# VarScan calls all sites with at least 5X coverage, a variant frequency above
# 10%, and 95% probability of being a SNP. Need varscan 2.3.5 version
java -jar VarScan.jar mpileup2snp mpileup --output-vcf --min-coverage 5 --strand-filter 0 --min-var-freq 0.1 --p-value 0.05 >SNPS.vcf
# VCFtools to filter raw SNPs and create a filtered vcf file (Final.recode.vcf)
# with SNPs that are present in every individual and that are not INDels
# can also work with --geno 0.99 flag but it needs vcftools 0.1.10 version
vcftools --vcf SNPS.vcf  --out Final --recode --non-ref-af 0.001 --remove-indels
# VCFtools again to filter for SNPs that are present at an average of 10X coverage 
vcftools --vcf Final.recode.vcf --out Final10X --recode --min-meanDP 10 
