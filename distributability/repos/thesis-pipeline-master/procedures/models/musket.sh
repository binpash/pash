# $MUSKET

# MUSKET (version 1.1) is a parallel multi-stage k-mer based error corrector
# Usage: musket [options] file [file1 ...]
# Basic Options:
# 	-k <int uint> (specify two paramters: k-mer size and estimated total number of k-mers for this k-mer size)
# 	   (e.g. estimated number of k-mers: 67108864, 134217728, 268435456 and 536870912)
# 	-o <str> (the single output file name)
# 	-omulti <str> (prefix of output file names, one input corresponding one output)
# 	-p <int> (number of threads [>=2], default=2)
# 	-zlib <int> (zlib-compressed output, default=0)
# 	-maxtrim <int> (maximal number of bases that can be trimmed, default=0)
# 	    (keeping the longest error-free substring of a read
# 	-inorder (keep sequences outputed in the same order with the input)
# 	-lowercase (write corrected bases in lowercase, default=0)
# Advanced:
# 	-maxbuff <int> (capacity of message buffer for each worker, default=1024)
# 	-multik <bool> (enable the use of multiple k-mer sizes, default=0)
# 	-maxerr <int> (maximal number of mutations in any region of length #k, default=4)
# 	-maxiter <int> (maximal number of correcting iterations per k-mer size, default=2)
# 	-minmulti <int> (minimum multiplicty for correct k-mers [only applicable when not using multiple k-mer sizes], default=0)
