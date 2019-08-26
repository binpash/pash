#! /home/users/cordier/.linuxbrew/bin/python3

if __name__ == "__main__":

	# Imports
	import sys, argparse
	# Library Import
	from Bio import SeqIO

	# 
	# Parse Arguments
	# 

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type = str, help = "Input FastQ File")
	parser.add_argument("-o", "--outputprefix", type = str, help = "Prefix to Output Fasta & Qual Files")
	argsDict = vars(parser.parse_args())

	fastq = argsDict["input"]
	prefix = argsDict["outputprefix"]

	# Assertions for Required Input
	assert (fastq is not None), "No FASTQ input provided!"

	#  If No Prefix, Use Same as FastQ
	if prefix is None: 
		prefix = ".".join(fastq.split(".")[0:-1])

	# 
	# Conversion
	#

	# Split FastQ into Fasta & Qual
	SeqIO.convert(fastq, "fastq", prefix + ".fasta", "fasta")
	print("Fasta file written to %s.fasta" % prefix)
	SeqIO.convert(fastq, "fastq", prefix + ".qual", "qual")
	print("Qual file written to %s.qual" % prefix)

else:

	pass
