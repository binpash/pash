#! /home/users/cordier/.linuxbrew/bin/python3

if __name__ == "__main__":

	# Imports
	import sys, argparse
	# Library Import
	from Bio import SeqIO
	from Bio.SeqIO.QualityIO import PairedFastaQualIterator

	# 
	# Parse Arguments
	# 

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--inputfasta", type = str, help = "Input Fasta File")
	parser.add_argument("-q", "--inputqual", type = str, help = "Input Qual File")
	parser.add_argument("-o", "--outputprefix", type = str, help = "Prefix to Output FastQ File")
	argsDict = vars(parser.parse_args())

	fasta = argsDict["inputfasta"]
	qual = argsDict["inputqual"]
	prefix = argsDict["outputprefix"]

	# Assertions for Required Input
	assert (fasta is not None), "No Fasta input provided!"
	assert (qual is not None), "No Qual input provided!"

	#  If No Prefix, Use Same as FASTQ
	if prefix is None: 
		prefix = ".".join(fasta.split(".")[0:-1])
		
	# 
	# Conversion
	#

	# Merge Fasta & Qual into FastQ
	records = PairedFastaQualIterator(open(fasta), open(qual))
	SeqIO.write(records, prefix + ".fastq", "fastq")

else:

	pass
