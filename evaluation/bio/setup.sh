mkdir -p input
mkdir -p output
cd input
mkdir -p bio3
mkdir -p bio2
wget ndr.md/data/bio/{R1.fastq.gz,R2.fastq.gz,ref.fa}
gunzip R1.fastq.gz
gunzip R2.fastq.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz 
gunzip human_g1k_v37.fasta.gz                                                             
cp R1.fastq bio2/Sample1.R1.fq
cp R2.fastq bio2/Sample1.R2.fq
git clone https://github.com/lh3/wgsim
cd wgsim/ && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm && cd -
apt-get install samtools bowtie2 vcftools sra-toolkit cutadapt
wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.zip
unzip 0.6.6.zip
mv TrimGalore-0.6.6 Trim/
wget https://phoenixnap.dl.sourceforge.net/project/bio-rainbow/rainbow_2.0.4.tar.gz
tar xf rainbow_2.0.4.tar.gz
mv rainbow_2.0.4 rainbow
cd rainbow && make -j && cd -
wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
tar xf fastx_toolkit-0.0.14.tar.bz2
wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz
tar xf libgtextutils-0.7.tar.gz
# patch because they suck
sed -i 's/return input_stream/return (bool)input_stream/g' ./libgtextutils-0.7/src/gtextutils/text_line_reader.cpp
cd libgtextutils-0.7 && ./configure && sudo make install && cd - 
sed -i 's/usage();/usage();break;/g' fastx_toolkit-0.0.14/src/fasta_formatter/fasta_formatter.cpp
cd fastx_toolkit-0.0.14  && ./configure && sudo make install && cd -
wget https://altushost-swe.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar xf bwa-0.7.17.tar.bz2
cd bwa-0.7.17 && make -j && cd -
wget https://deac-fra.dl.sourceforge.net/project/varscan/VarScan.v2.3.5.jar -O bio2/VarScan.jar
#### Download sam files ###
wget https://raw.githubusercontent.com/samtools/samtools/develop/examples/toy.sam
#wget https://raw.githubusercontent.com/samtools/samtools/develop/examples/ex1.sam.gz
#gunzip ex1.sam.gz
#sam-dump --aligned-region 20:1-16444167 --output-file SRR1976040_chr20.sam SRR1976040
# too slow to download
#sam-dump --aligned-region chr20 --output-file SRR1976036_chr20.sam SRR1976036

## bio3.sh
cd bio3
wget www.bioinformaticsworkbook.org/Appendix/GNUparallel/fastqfiles.tar.gz
tar -zxvf fastqfiles.tar.gz
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
cd -
