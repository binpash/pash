#!/bin/bash

#  Identify the top 10 reasons why  genome assemblies don't make it into GenBank
# -- NIH's  genetic sequence database,  an annotated collection of  all publicly
# available DNA sequences
# http://thegenomefactory.blogspot.com/2019/09/25-reasons-assemblies-dont-make-it-into.html

# Require: csvkit
# Data: http://ndr.md/data/bio/genbank.txt

IN=./input/genbank.txt
OUT=./output/out.txt

cat $IN |
  csvcut -t -K 1 -c 'excluded_from_refseq' |
  tail -n +2 | tr ";" "\n" |
  sed -e 's/^ //' -e 's/ $//' |
  grep -v '""' |
  sort |
  uniq -c |
  sort -nr |
  head -n 10 |
  nl > $OUT

# More bio pipelines for 
# # Strains with Complete Genome
# cat assembly_summary.tsv \
#     | csvtk grep -t -f assembly_level -i -p "Complete Genome"  \
#     | wc -l
# 
# # Most sequenced species with Complete Genome
# cat assembly_summary.tsv \
#     | csvtk grep -t -f assembly_level -i -p "Complete Genome"  \
#     | csvtk cut -t -f organism_name \
#     | cut -d ' ' -f 1,2 \
#     | csvtk freq -t -n -r | head -n 20 | csvtk pretty -t
# 
# # Number of species, by organism name
# 
# # Filter by species (organism_name)
# cat assembly_summary.tsv \
#     | csvtk grep -t -f organism_name -i -r -p "Mycobacterium tuberculosis" \
#     | csvtk grep -t -f assembly_level -i -p "Complete Genome"  \
#     > mt.tsv
# 
# # Filter (complete genome) by species_taxid
# cat assembly_summary.tsv \
#     | csvtk grep -t -f species_taxid -p 239935,1280 \
#     | csvtk grep -t -f assembly_level -i -p "Complete Genome" \
#     > bytaxid.tsv
# 
# # Download genome sequence and annotation files
# cat mt.tsv | csvtk cut -t -f ftp_path | sed 1d \
#     | rush -v prefix='{}/{%}' \
#         ' \
#             wget -c {prefix}_genomic.fna.gz; \
#             wget -c {prefix}_genomic.gbff.gz; \
#             wget -c {prefix}_genomic.gff.gz; \
#             wget -c {prefix}_cds_from_genomic.fna.gz \
#             wget -c {prefix}_protein.faa.gz; \
#         ' \
#         -j 10 -c -C download.rush
# 
# #Get GenBank assembly summary file
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
# 
# #Get all lines that have "Mycobacter", if 12th field is "Complete Genome", print the 20th field (url to file).
# #But the actual filename ends _genomic.fna.gz so include that too..
# grep Mycobacter assembly_summary_genbank.txt \
#     | awk 'BEGIN{FS="\t"}{if($12=="Complete Genome"){print $20}}' \
#     | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.fna.gz"}' \
#     > urls.txt
# 
# #Now you can go through your urls file
# IFS=$'\n'; for NEXT in $(cat urls.txt); do wget "$NEXT"; done
