
# $DECGPU
# DecGPU is a distributed error correction on GPUs using CUDA and MPI (CPU-based)
# Usage:
# 	decgpu directory [options] [-fasta | -fastq ] infile1 [infile2]

# 	directory	working directory
# Options:
# 	-fasta infile1 [infile2] (input reads file in FASTA format)
# 	-fastq infile1 [infile2] (input reads file in FASTQ format)
# 	-k <integer> (the kmer size (odd number between 0 and 32), default value: 21 )
# 	-paired (indicating all the input reads are paired-end internally, and the output are also paired)
# 	-numthreads <integer> (the number of OpenMP threads per process, default value: 2 )
# 	-minmulti <integer> (the minimum multiplicty cutoff, default value: 6 )
# 	-minvotes <integer> (the mininum votes per nucleotide per position, default value: 3 )
# 	-numsearch <integer> (the number of searches, default value: 1 )
# 	-maxtrim <integer> (the maximum number of bases that can be trimmed, default value: 4 )
# 	-est_bf_size (estimate the bloom filter size (recommended for small datasets, use the MAX size by default)
# 	-version (print out the version)