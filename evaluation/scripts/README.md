## Summary of Results

This table is for keeping track of `Dish` results, and blocking aspects when `Dish` can't execute something distributed. All input is provided in [./input](./input) (possibly requires executing `./input/fetch.sh`). 

| Pipeline                         | 1x                                  | 100x    | Notes                                         |
| -------------------------------- | ----------------------------------- | ------- | --------------------------------------------- |
| [genome-diff.sh](genome-diff.sh) | 101.98s, 4.88s, 114% — _1:32.99_    |         |                                               |
| [diff.sh](diff.sh)               | 452.90s, 85.40s, 234% — _3:49.10_   |         |                                               |
| [set-diff.sh](set-diff.sh)       | 305.81s, 3.35s, 204% — _2:31.26_    |         |                                               |
| [genquality.sh](genquality.sh)   | 7.03s, 2.30s, 148% — _6.305_        |         |                                               |
| [search.sh](search.sh)           | 4825.71s, 5.93s, 99% — _1:20:33.88_ |         |                                               |
| [bigrams.sh](bigrams.sh)         | 738.01s, 34.73s, 107% — _11:59.99_  |         |                                               |
| [trigrams.sh](trigrams.sh)       | 971.30s, 60.17s, 111% — _15:29.22_  |         |                                               |
| [page-count.sh](page-count.sh)   | 0.21s, 0.03s, 107% — _0.228_        |         |                                               |
| [tailprogs.sh](tailprogs.sh)     | 1.53s, 0.34s, 103% — _1.816_        |         |                                               |
| [spell.sh](spell.sh)             | 0.74s, 0.07s, 176% — _0.462_        |         |                                               |
| [symtab-sha.sh](symtab-sha.sh)   | 0.43s, 0.21s, 362% — _0.177_        |         |                                               |
| [topn.sh](topn.sh)               | 422.06s, 12.94s, 105% — _6:53.29_   |         |                                               |
| [wc.sh](wc.sh)                   | 14.07s, 5.54s, 137% — _14.289_      |         |                                               |
| [wf.sh](wf.sh)                   | 129.66s, 10.00s, 106% — _2:10.59_   |         |                                               |
| [longest-man.sh](longest-man.sh) | 0.09s, 0.14s, 150% — _0.151_        |         |                                               |
|                                  |                                     |         |                                               |
| [compile.sh](compile.sh)         | 18.02s, 0.50s, 100% — _18.489_      |         | Does not make sense to parallelize            |
| [sieve.sh](sieve.sh)             | 12.41s, 8.78s, 68% — _30.816_       |         | Not a pipeline                                |
| [maximal.sh](maximal.sh)         |                                     |         | Only for compatibility testing                |
| [merge-wc.sh](merge-wc.sh)       |                                     |         | Tool                                          |
| [old_bigrams.sh](old_bigrams.sh) |                                     |         | Do we need it?                                |

All scripts write to `./output/out.txt` (but could be made to read from and write to `/dev/shm/out.txt`).

To run and time all scripts, execute:

```sh
./input/fetch.sh
for f in *.sh; do
  time ./$f
done
```

## Pipeline Details

This directory contains a few example pipelines:

* [genome-diff.sh](genome-diff.sh): Find differences between two genome sequences---a a paired Illumina sequencing read  (FASTQ files)  and an  assembled  reference genome  from GenBank  (e.g., Pasteurella multocida).
* [compile.sh](compile.sh): Find markdown files  in the current directory tree, compile  them to HTML, and serve them over the network.
* [diff.sh](diff.sh): Compare two shuffled-and-sorted streams, element by element.
* [set-diff.sh](set-diff.sh): Show the set-difference between two streams (i.e., elements in the first that are not in the second).
* [genquality.sh](genquality.sh): Identify the primary reasons why GenBank rejects genome assemblies.
* [search.sh](search.sh): Substring search of a complex string, using backtracking (equivalent of "Hello, World!" Program).
* [maximal.sh](maximal.sh): A script for testing parsing compatibility.
* [merge-wc.sh](merge-wc.sh): (**should move to a `./lib/wrappers`?**)
* [bigrams.sh](bigrams.sh): Identify all bigrams in an input text.
* [bigrams.sh](bigrams.sh): Identify all trigrams in an input text.
* [old_bigrams.sh](old_bigrams.sh): Identify all bigrams in an input text (using interm. files).
* [page-count.sh](page-count.sh): Determine how many pages are in a folder of OpenOffice documents.
* [tailsprogs.sh](tailsprogs.sh): Find the 25 longest programs/shell scripts in the system.
* [sieve.sh](sieve.sh): Sieve of Eratosthenes for finding all prime numbers up to a limit.
* [spell.sh](spell.sh): Spell-check one or more man-pages.
* [symtab-sha.sh](symtab-sha.sh): Extract and hash an executable's symbol, useful for building and signing an SGX enclave
* [topn.sh](topn.sh): Identify top 10 most frequent terms in an input text stream.
* [wc.sh](wc.sh): Count all words in a stream.
* [wf.sh](wf.sh): Report on the word-frequencies of distinct terms in an input text stream.
* [longest-man.sh](longest-man.sh): Reports on the longest manual pages found in the system.

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
cd ~
BIO_TOOLS=~/biotools

# These should be added to every script
export PATH="$PATH:$BIO_TOOLS/bcftools-1.9"
export PATH="$PATH:$BIO_TOOLS/samtools-1.9"
export PATH="$PATH:$BIO_TOOLS/htslib-1.9"
export PATH="$PATH:$BIO_TOOLS/minimap2-2.17_x64-linux"

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

cd ..
curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf - 
cd ..


pip3 install csvkit

# page-count.sh
sudo apt-get install libimage-exiftool-perl
sudo apt-get install bc
```
