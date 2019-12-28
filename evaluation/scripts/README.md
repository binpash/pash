## Summary of Results

This table is for keeping track of `Dish` results, and blocking aspects when `Dish` can't execute something.

| Pipeline              | 1x     | 100x    | Notes    |
| --------------------- | ------ | ------- | -------- |
| [genome-diff.sh]      |        |         |          |         
| [compile.sh]          |        |         |          |     
| [diff.sh]             |        |         |          |  
| [set-diff.sh]         |        |         |          |      
| [genquality.sh]       |        |         |          |        
| [grep.sh]             |        |         |          |  
| [maximal.sh]          |        |         |          |     
| [merge-wc.sh]         |        |         |          |      
| [minimal.sh]          |        |         |          |     
| [bigrams.sh]          |        |         |          |     
| [old_bigrams.sh]      |        |         |          |         
| [page-count.sh]       |        |         |          |        
| [tailprogs.sh]        |        |         |          |       
| [sieve.sh]            |        |         |          |   
| [spell.sh]            |        |         |          |   
| [symtab-sha.sh]       |        |         |          |        
| [topn.sh]             |        |         |          |  
| [wc.sh]               |        |         |          |
| [wf.sh]               |        |         |          |


## Pipeline Details

This directory contains a few example pipelines:

* [genome-diff.sh]: Find differences between two genome sequences---a a paired Illumina sequencing read  (FASTQ files)  and an  assembled  reference genome  from GenBank  (e.g., Pasteurella multocida).
* [compile.sh]: Find markdown files  in the current directory tree, compile  them to HTML, and serve them over the network.
* [diff.sh]: Compare two shuffled-and-sorted streams, element by element.
* [set-diff.sh]: Show the set-difference between two streams (i.e., elements in the first that are not in the second).
* [genquality.sh]: Identify the primary reasons why GenBank rejects genome assemblies.
* [grep.sh]: Substring search of a complex string, using backtracking.
* [maximal.sh]: A script for testing parsing compatibility.
* [merge-wc.sh]: (**should move to a `./lib/wrappers`?**)
* [minimal.sh]: A small pipeline that contains only `grep` and `tr`
* [bigrams.sh]: Identify all bigrams in an input text.
* [old_bigrams.sh]: Identify all bigrams in an input text (using interm. files).
* [page-count.sh]: Determine how many pages are in a folder of OpenOffice documents.
* [tailsprogs.sh]: Find the 25 longest programs/shell scripts in the system.
* [sieve.sh]: Sieve of Eratosthenes for finding all prime numbers up to a limit.
* [spell.sh]: Spell-check one or more man-pages.
* [symtab-sha.sh]: Extract and hash an executable's symbol, useful for building and signing an SGX enclave
* [topn.sh]: Identify top 1000 most frequent terms in an input text stream.
* [wc.sh]: Count all words in a stream.
* [wf.sh]: Report on the word-frequencies of distinct terms in an input text stream.

These files are candidates for removal:

* grep2.sh 
* minimal2.sh
* minimal3.sh
* minimal4.sh
* minimal5.sh
* pretty_print_json.sh* -- move to a `/tools` dir or simply pipe to `jq .`?
* call_distrib_planner_example.sh -- move to a `/tools` dir?

## Setup Details

```sh
# Install dependencies for some tools
BIO_TOOLS="~/biotools/"
mkdir $BIO_TOOLS

sudo apt-get update
sudo apt-get install gcc
sudo apt-get install make
sudo apt-get install libbz2-dev
sudo apt-get install zlib1g-dev
sudo apt-get install libncurses5-dev 
sudo apt-get install libncursesw5-dev
sudo apt-get install liblzma-dev

cd $BIO_TOOLS
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
make

cd ..
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make

cd ..
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
make

export PATH="$PATH:/usr/bin/bcftools-1.9"
export PATH="$PATH:/usr/bin/samtools-1.9"
export PATH="$PATH:/usr/bin/htslib-1.9"
source ~/.profile

pip3 install csvkit

# page-count.sh
sudo apt install libimage-exiftool-perl
sudo apt install bc
```
