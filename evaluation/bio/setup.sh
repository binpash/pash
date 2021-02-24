cd $PASH_TOP/evaluation/scripts/input
wget ndr.md/data/bio/{R1.fastq.gz,R2.fastq.gz,ref.fa}
gunzip R1.fastq.gz
gunzip R2.fastq.gz
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz 
#gunzip human_g1k_v37.fasta.gz                                                             
cp R1.fastq Sample1.R1.fq
cp R2.fastq Sample1.R2.fq
git clone https://github.com/lh3/wgsim
cd wgsim/ && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm && cd -
#wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
#tar xf samtools-1.11.tar.bz2
apt-get install samtools bowtie2  vcftools
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
wget https://deac-fra.dl.sourceforge.net/project/varscan/VarScan.v2.3.5.jar -O VarScan.jar
