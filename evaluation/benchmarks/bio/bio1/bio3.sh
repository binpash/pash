# **Create the bowtie2 alignment database for the Arabidopsis genome**
# https://bioinformaticsworkbook.org/Appendix/GNUparallel/GNU_parallel_examples.html#gsc.tab=0
cd $PASH_TOP/evaluation/bio/input/bio3
bowtie2-build TAIR10_chr_all.fas tair
#theirs                                                                          
time parallel -j2 "bowtie2 --threads 4 -x tair -k1 -q -1 {1} -2 {2} -S {1/.}.sam  >& {1/.}.log" ::: fastqfiles/*_1.fastq.gz :::+ fastqfiles/*_2.fastq.gz
#ours                                                                            
paste  <(find . -name "*_1.fastq.gz") <(find . -name "*_2.fastq.gz") | xargs -n \
2 sh -c 'bowtie2 --threads 4 -x tair -k1 -q -1 "$1" -2 "$2" -S fifth_R1.sam' argv0

