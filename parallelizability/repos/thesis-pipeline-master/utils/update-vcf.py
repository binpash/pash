#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
A Quick & Dirty Script to Update Chromosome Names in a VCF (e.g. From 'chr1' to '1')
"""

import re       as Regex
import sys      as System

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

if __name__ == "__main__":
	# Generate Replacement Dict
	updateDict      = dict()
	chroms          = [i + 1 for i in range(22)]
	chroms.extend(["X", "Y", "M", "Un"])
	for chrID in chroms:
		chrID = str(chrID)
		updateDict["/chr" + chrID + "/"] = "/chr" + chrID + "/"
		updateDict[".chr" + chrID + "."] = ".chr" + chrID + "."
		updateDict["chr" + chrID + "\t"] = chrID + "\t"
		updateDict["=chr" + chrID + ","] = "=" + chrID + ","
		updateDict["=chr" + chrID + "_"] = "=" + chrID + ","

	files = System.argv[1:]
	i = 1
	for file in files:
		name = file.replace(".vcf", "")
		with open(file, "r") as VCF:
			with open(name + "_modified.vcf", "w") as updateVCF:
				for line in VCF:
					for key, value in updateDict.items():
						line = line.replace(key, value)
					updateVCF.write(line)
		print(str(i), "of", str(len(files)))
		i += 1
	print("DONE")
else:
	pass

