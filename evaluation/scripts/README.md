## Summary of Results

This table is for keeping track of `Dish` results, and blocking aspects when `Dish` can't execute something distributed. All input is provided in [./input](./input) (possibly requires executing `./input/fetch.sh`). 

| Pipeline              | 1x                                                  | 100x    | Notes                                         |
| --------------------- | --------------------------------------------------- | ------- | --------------------------------------------- |
| [genome-diff.sh]      | 101.98s user 4.88s system 114% cpu 1:32.99 total    |         |                                               |
| [diff.sh]             | 452.90s user 85.40s system 234% cpu 3:49.10 total   |         |                                               |
| [set-diff.sh]         | 305.81s user 3.35s system 204% cpu 2:31.26 total    |         |                                               |
| [genquality.sh]       | 7.03s user 2.30s system 148% cpu 6.305 total        |         |                                               |
| [search.sh]           | 4825.71s user 5.93s system 99% cpu 1:20:33.88 total |         |                                               |
| [bigrams.sh]          | 738.01s user 34.73s system 107% cpu 11:59.99 total  |         |                                               |
| [trigrams.sh]         | 971.30s user 60.17s system 111% cpu 15:29.22 total  |         |                                               |
| [page-count.sh]       | 0.21s user 0.03s system 107% cpu 0.228 total        |         |                                               |
| [tailprogs.sh]        | 1.53s user 0.34s system 103% cpu 1.816 total        |         |                                               |
| [spell.sh]            | 0.74s user 0.07s system 176% cpu 0.462 total        |         |                                               |
| [symtab-sha.sh]       | 0.43s user 0.21s system 362% cpu 0.177 total        |         |                                               |
| [topn.sh]             | 422.06s user 12.94s system 105% cpu 6:53.29 total   |         |                                               |
| [wc.sh]               | 14.07s user 5.54s system 137% cpu 14.289 total      |         |                                               |
| [wf.sh]               | 129.66s user 10.00s system 106% cpu 2:10.59 total   |         |                                               |
| [longest-man.sh]      | 0.09s user 0.14s system 150% cpu 0.151 total        |         |                                               |
|                       |                                                     |         |                                               |
| [compile.sh]          | 18.02s user 0.50s system 100% cpu 18.489 total      |         | Does not make sense to parallelize            |
| [sieve.sh]            | 12.41s user 8.78s system 68% cpu 30.816 total       |         | Not a pipeline                                |
| [maximal.sh]          |                                                     |         | Only for compatibility testing                |
| [merge-wc.sh]         |                                                     |         | Tool                                          |
| [old_bigrams.sh]      |                                                     |         | Do we need it?                                |

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

* [genome-diff.sh]: Find differences between two genome sequences---a a paired Illumina sequencing read  (FASTQ files)  and an  assembled  reference genome  from GenBank  (e.g., Pasteurella multocida).
* [compile.sh]: Find markdown files  in the current directory tree, compile  them to HTML, and serve them over the network.
* [diff.sh]: Compare two shuffled-and-sorted streams, element by element.
* [set-diff.sh]: Show the set-difference between two streams (i.e., elements in the first that are not in the second).
* [genquality.sh]: Identify the primary reasons why GenBank rejects genome assemblies.
* [search.sh]: Substring search of a complex string, using backtracking (equivalent of "Hello, World!" Program).
* [maximal.sh]: A script for testing parsing compatibility.
* [merge-wc.sh]: (**should move to a `./lib/wrappers`?**)
* [bigrams.sh]: Identify all bigrams in an input text.
* [bigrams.sh]: Identify all trigrams in an input text.
* [old_bigrams.sh]: Identify all bigrams in an input text (using interm. files).
* [page-count.sh]: Determine how many pages are in a folder of OpenOffice documents.
* [tailsprogs.sh]: Find the 25 longest programs/shell scripts in the system.
* [sieve.sh]: Sieve of Eratosthenes for finding all prime numbers up to a limit.
* [spell.sh]: Spell-check one or more man-pages.
* [symtab-sha.sh]: Extract and hash an executable's symbol, useful for building and signing an SGX enclave
* [topn.sh]: Identify top 10 most frequent terms in an input text stream.
* [wc.sh]: Count all words in a stream.
* [wf.sh]: Report on the word-frequencies of distinct terms in an input text stream.
* [longest-man.sh]: Reports on the longest manual pages found in the system.

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
