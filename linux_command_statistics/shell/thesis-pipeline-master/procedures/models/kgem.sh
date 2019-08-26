# java -jar $KGEM
# Error: Missing argument ReadsFile
# Error: Missing argument k
# kGEM version 0.3.1: Local Reconstruction for Mixed Viral Populations.
# Usage: kGEM [options] ReadsFile k

#   ReadsFile
#         Fasta or Sam file containing aligned sequence data.
#   k
#         Number of initial seeds for clustering. May be a single positive
# 	integer value or may be a range to try for model selection (e.g. 5:10).
#   -t <value> | --threshold <value>
#         Min Hamming distance between seeds for init stage.
#   -f <value> | --scoring-func <value>
#         The scoring method to use in model selection. Only used when a
# 	range is specified. (Default: 'AICc').
# 	Accepted Values:

# 		AIC: Akaike Information Criteria
# 		AICc: Corrected Akaike Information Criteria
# 		BIC: Bayesian Information Criteria
# 		BICMAP: Bayesian Information Criteria using MAP estimate
#   -g <value> | --consensus-file <value>
#         Optional Fasta File containing initial seeds.
#   -o <value> | --output-dir <value>
#         Directory to output results.
#   -h | --help
#         Prints this help text.
#   -v | --version
#         Prints kGEM version.