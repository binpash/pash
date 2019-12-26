#!/bin/bash

# also, with the same dataset:
# https://www.biostars.org/p/255212/


pip3 install csvkit

csvcut -t -K 1 -c 'excluded_from_refseq' <(curl -s https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt) |
  tail -n +2 | tr ";" "\n" |
  sed -e 's/^ //' -e 's/ $//' |
  grep -v '""' |
  sort |
  uniq -c |
  sort -nr 
