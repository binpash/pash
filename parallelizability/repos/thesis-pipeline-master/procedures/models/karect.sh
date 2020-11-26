# $KARECT -correct
# -Run "karect -correct [options_list]" for the error correction tool.
# -Available options: (i=integer, d=decimal, f=file, s=directory)
# -Essential options:
#   "-inputfile=f":      Specify an input fasta/fastq file. This option can be repeated for multiple files
#   "-celltype=[haploid|diploid]": Specify the cell type. Use "haploid" for bacteria and viruses.
#   "-matchtype=[edit|hamming|insdel]": Specify the matching type. "hamming" allows substitution errors only.
#                          "edit" allows insertions, deletion, and substitutions with equal costs.
#                          "insdel" is the same as "edit", but the cost of substitutions is doubled.
#                          Use "hamming" for Illumina datasets, and "edit" for 454 datasets.
# -Basic options:
#   "-inputdir=s":       Specify the input directory. Ignored if input file paths are complete [Default=.].
#   "-resultdir=s":      Specify a directory to save result file(s) [Default=.].
#   "-resultprefix=s":   Specify a prefix string to the result file(s) [Default=karect_].
#   "-tempdir=s":        Specify a directory to save temporary output files [Default=.].
#   "-threads=i":        Specify the number of threads [Default=16].
#   "-memory=d":         Specify an upper bound on the memory that can be used in gigabytes [Default=10240.0].
# -Advanced options:
#   "-aggressive=d":     Specify the aggressiveness towards error correction [Default=0.42].
#   "-numstages=i":      Specify the number of stages (1 or 2) [Default=1].
#   "-minoverlap=i":     Specify the minimum overlap size [Default=35].
#   "-minoverlapper=d":  Specify the minimum overlap percentage [Default=0.20].
#   "-minreadweigth=d":  Specify the minimum read weight [Default=0].
#   "-errorrate=d":      Specify the first stage maximum allowed error rate [Default=0.25].
#   "-errorratesec=d":   Specify the second stage maximum allowed error rate [Default=0.25].
#   "-reserveval=d":     Specify the minimum reservation value [Default=100.0].
#   "-estcov=[yes|no]":  Estimate coverage and use it to adjust the minimum reservation value [Default=yes].
#   "-usequal=[yes|no]": Use quality values of candidate reads [Default=yes].
#   "-higherror":        Work in high error rate mode (error rate = 0.50).
#   "-trimfact=d":       Specify the trimming factor [Default=2.5].
#   "-reserveper=d":     Specify the minimum reservation percentage [Default=1.0].
#   "-kmer=i":           Specify the minimum kmer size (will increase according to '-kmerfactor') [Default=9].
#   "-trim=[yes|no]":    Allow/Disallow trimming. Do not allow trimming if evaluating results afterwards, or
#                          if you will pass results to an assembler which expects fixed read sizes [Default=no].
# -More advanced options:
#   "-maxlenmatches=i":  Specify the maximum number of expected alignment computations (millions) [Default=2000].
#   "-maxkmerslots=i":   Specify the maximum number kmer slots to be used [Default=100,000].
#   "-kmerfactor=i":     Specify the factor f such that 4^kmersize > total_num_kmers/f [Default=1000].
#   "-maxkmerres=i":     Specify the maximum number kmer results to be used [Default=30].
#   "-kmererrors=i":     Specify the maximum allowed kmer errors (0,1,2) [Default=2].
#   "-readsperstep=i":   Specify the maximum number of processed reads per step [Default=1000].
#   "-fbs=i":            Specify file block size in megabytes [Default=10].
#   "-cbs=i":            Specify cache block size in megabytes [Default=128].
# -Example:
#         ./karect -correct -inputdir=/sra_data -inputfile=SRR001666_1.fasta -inputfile=SRR001666_2.fasta
#               -resultprefix=karect_r1_ -celltype=haploid -matchtype=hamming -errorrate=0.25 -threads=12