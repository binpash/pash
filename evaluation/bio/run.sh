# Can we run this or something like this?
# https://github.com/marcelm/cutadapt/issues/157
# https://www.biostars.org/p/123237/
# Some fastq files 
# wget ndr.md/data/bio/{R1.fastq.gz,R2.fastq.gz,ref.fa}

PW=$PASH_TOP/evaluation/scripts/input
#remove adapter
remove_adapter()
(
  find . -name "$PW/*.fastq" |  sort | uniq | xargs -I {} cutadapt -a AGATCGGAAGAGCACAC {} >  /dev/null
)

# convert fastq to fasta format
# It recognizes the extension .fasta and it converts the input to tha format
convert_to_fasta()
(
  find . -name "$PW/*.fastq" | xargs -I {} cutadapt -o {}.fasta.gz {}
)
# remove more than once adapter, run the tool twice
# We could create random adapter inputs for several fastq
remove_adapter_twice()
(
  cutadapt -g ^TTAAGGCC -g ^AAGCTTA $PW/R1.fastq | cutadapt -a TACGGACT - > /dev/null
  cutadapt -g ^TTAAAACC -g ^AAGCTTA $PW/R1.fastq | cutadapt -a TACGAACT - > /dev/null
)

# trim primers
trim_primers()
(
  find . -name "$PW/*.fastq" | xargs -I {}  cutadapt -a TCCTCCGCTTATTGATAGC -o ${i}\_trimmed.fastq {}; 
)

# convert sam to bam
# Need INPUT HERE
sam_to_bam()
(
  
  for i in "${NAMES[@]}"
  do
    samtools view -bS -q15 $i.sam > $i.bam
    samtools sort $i.bam $i
  done 
)

# Here are sample steps to generate a single paired read from hg19:
# https://www.biostars.org/p/150010/
pipeline()
(
  PW=$PASH_TOP/evaluation/scripts/input
  cd $PW
  echo $PW
  #download hg19 reference genome, e.g.
  #wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
  #wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
  #gunzip human_g1k_v37.fasta.gz

  #filter out a single chromosome and index it, e.g.
  samtools faidx $PW/human_g1k_v37.fasta 20 > $PW/human_g1k_v37_chr20.fasta
  bowtie2-build $PW/human_g1k_v37_chr20.fasta $PW/homo_chr20
  #simulate a single read sample, e.g. here is for a single (-N 1) paired read:
  wgsim/wgsim -N 1 $PW/human_g1k_v37_chr20.fasta $PW/single.read1.fq $PW/single.read2.fq > $PW/wgsim.out
  #generate the sam, e.g.
  bowtie2 -x $PW/homo_chr20 -1 $PW/single.read1.fq -2 $PW/single.read2.fq -S $PW/single_pair.sam
  #generate a bam
  samtools view -b -S -o $PW/single_pair.bam $PW/single_pair.sam 
  #sort and index it
  samtools sort $PW/single_pair.bam -o $PW/single_pair.sorted.bam
  # this seems to not affect the file, but in other cases, its indeed needed
  samtools index $PW/single_pair.sorted.bam 
)

# https://dfzljdn9uc3pi.cloudfront.net/2013/203/1/Supplement_S2.pdf
# #Script to automatically process X number samples, produce a reference
# assembly, map reads back to the assembly, and call SNPs 
# https://zenodo.org/record/940733


