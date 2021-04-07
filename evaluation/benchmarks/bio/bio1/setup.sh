if [[ $1 == "-c" ]]; then
  rm -rf input
  rm -rf output
  exit
fi

mkdir -p input
mkdir -p output
cd input
if [[ ! -f R1.fastq ]]; then
  wget ndr.md/data/bio/{R1.fastq.gz,R2.fastq.gz,ref.fa}

  gunzip R1.fastq.gz
  gunzip R2.fastq.gz
fi
if [[ ! -f human_g1k_v37.fasta ]]; then
  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz 
  gunzip human_g1k_v37.fasta.gz                                                             
  cp R1.fastq Sample1.R1.fq
  cp R2.fastq Sample1.R2.fq
fi

if [[ ! -d wgsim ]]; then
  git clone https://github.com/lh3/wgsim
  cd wgsim/ && gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm && cd -
fi
apt-get install samtools bowtie2 vcftools sra-toolkit cutadapt zlib1g-dev default-jre
if [[ ! -d Trim ]]; then
  wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.zip
  unzip 0.6.6.zip
  mv TrimGalore-0.6.6 Trim/
fi
if [[ ! -d rainbow ]]; then
  wget https://phoenixnap.dl.sourceforge.net/project/bio-rainbow/rainbow_2.0.4.tar.gz
  tar xf rainbow_2.0.4.tar.gz
  mv rainbow_2.0.4 rainbow
  cd rainbow && make -j && cd -
fi
if [[ ! -d fastx_toolkit-0.0.14 ]]; then
  wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
  tar xf fastx_toolkit-0.0.14.tar.bz2
  sed -i 's/usage();/usage();break;/g' fastx_toolkit-0.0.14/src/fasta_formatter/fasta_formatter.cpp
  cd fastx_toolkit-0.0.14  && ./configure && sudo make install && cd -
fi
if [[ -d libgtextutils-0.7 ]]; then
  wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz
  tar xf libgtextutils-0.7.tar.gz
  # patch because they suck
  sed -i 's/return input_stream/return (bool)input_stream/g' ./libgtextutils-0.7/src/gtextutils/text_line_reader.cpp
  cd libgtextutils-0.7 && ./configure && sudo make install && cd - 
fi 
if [[ ! -d bwa-0.7.17 ]]; then
  wget https://altushost-swe.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
  tar xf bwa-0.7.17.tar.bz2
  cd bwa-0.7.17 && make -j && cd -
fi

if [[ ! -f VarScan.jar ]]; then
  wget https://deac-fra.dl.sourceforge.net/project/varscan/VarScan.v2.3.5.jar -O VarScan.jar
fi

#### Download sam files ###
if [[ ! -f toy.sam ]]; then
  wget https://raw.githubusercontent.com/samtools/samtools/develop/examples/toy.sam
  wget https://raw.githubusercontent.com/samtools/samtools/develop/examples/ex1.sam.gz
  gunzip ex1.sam.gz
  sam-dump --aligned-region 20:1-16444167 --output-file SRR1976040_chr20.sam SRR1976040
  # too slow to download
  sam-dump --aligned-region chr20 --output-file SRR1976036_chr20.sam SRR1976036
fi

## bio3.sh
#wget www.bioinformaticsworkbook.org/Appendix/GNUparallel/fastqfiles.tar.gz
#tar -zxvf fastqfiles.tar.gz
#wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
#cd -
