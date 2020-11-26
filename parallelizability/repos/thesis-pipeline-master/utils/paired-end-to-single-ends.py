#! /home/users/cordier/.linuxbrew/bin/python3

def lineGenerator (filehandle, type = "r"):
    with open(filehandle, type) as file:
        for i, line in enumerate(file):
            yield line, i

if __name__ == "__main__":

    # Imports
    import argparse

    # 
    # Parse Arguments
    # 

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "Input Paired-End FastQ File")
    parser.add_argument("-o", "--outputprefix", type = str, help = "Prefix to Output Single End FastQ Files (i.e. prefix.1.fastq and prefix.2.fastq")
    argsDict = vars(parser.parse_args())

    fastq = argsDict["input"]
    prefix = argsDict["outputprefix"]

    # Assertions for Required Input
    assert (fastq is not None), "No FastQ input provided!"

    #  If No Prefix, Use Same as FASTQ
    if prefix is None: 
        prefix = ".".join(fastq.split(".")[0:-1])

    # File Handles
    handleA = prefix + ".1.fastq"
    handleB = prefix + ".2.fastq"

    # 
    # Conversion
    #

    # Open Singles, Open Paired, Iterate, & Write to Singles
    with open(handleA, "w") as singleA, open(handleB, "w") as singleB:
        lines = lineGenerator(fastq)
        for line, i in lines:
            if (i % 8 < 4):
                singleA.write(line)
            else:
                singleB.write(line)

else:

    pass
