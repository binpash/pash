#! /home/users/cordier/.linuxbrew/bin/python3

def interleave(single1, single2, strict = False):
    # import itertools
    for singleA, singleB in zip(single1, single2):
        if strict:
            assert singleA.id[0:-2] == singleB.id[0:-2], "Single-end IDs don't match (%s - %s)" % (singleA.id, singleB.id)
        singleA.id += "/1"
        singleB.id += "/2"
        yield singleA
        yield singleB

if __name__ == "__main__":

    # Imports
    import argparse
    # Library Imports
    from Bio import SeqIO

    # 
    # Parse Arguments
    # 
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", "--input1", type = str, help = "First Single-End FastQ File")
    parser.add_argument("-2", "--input2", type = str, help = "Second Single-End FastQ File")
    parser.add_argument("-o", "--outputprefix", type = str, help = "Prefix to Output Paired End FastQ File (i.e. prefix.fastq")
    parser.add_argument("-s", "--strict", action = "store_true", help = "Require all single-end IDs match their mate")
    argsDict = vars(parser.parse_args())

    handleA = argsDict["input1"]
    handleB = argsDict["input2"]
    prefix = argsDict["outputprefix"]

    # Assertions for Required Input
    assert (handleA is not None), "No first single-end FastQ input provided!"
    assert (handleB is not None), "No second single-end FastQ input provided!"

    # If No Prefix, Use Same as FASTQ
    if prefix is None: 
        if ".".join(singleA.split(".")[0:-1]) == ".".join(singleB.split(".")[0:-1]):
            prefix = ".".join(singleA.split(".")[0:-1])
            print("Output prefix: %s" % prefix)
    
    # If Single End Prefixes Are Not The Same
    assert (prefix is not None), "No output prefix provided (& input prefixes differ – this tool doesn't assume an output destination)!"

    # 
    # Conversion
    #

    # File Handles
    pairHandle = prefix + ".fastq"

    # Interleave Files & Write Output
    with open(pairHandle, "w") as paired:
        with open(handleA, "r") as singleA, open(handleB, "r") as singleB:
            recordsA = SeqIO.parse(singleA, "fastq")
            recordsB = SeqIO.parse(singleB, "fastq")
            SeqIO.write(interleave(recordsA, recordsB), paired, "fastq")

else:

    pass
