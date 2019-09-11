#! /home/users/cordier/.linuxbrew/bin/python3

#
# Note: This Script Has Only Been Validated on FastQ
#

def chunks (iterator, size):
    """
    Split File into Chunks for Processing - Not Yet Implemented
    """
    reads = True
    while reads:
        reads = []
        while len(chunk) < size:
            try:
                reads = iterator.next()
            except StopIteration:
                reads = None
            if reads is None:
                break
            chunk.append(reads)
        if chunk:
            yield chunk

if __name__ == "__main__":

    # Imports
    import sys, argparse
    # Library Import
    from Bio import SeqIO
    from Bio.SeqIO.QualityIO import PairedFastaQualIterator

    # Accepted Formats
    acceptedFormats = ["fasta", "fastq", "qual", "fa", "fq", "sam"]

    # 
    # Parse Arguments
    # 

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "Input File")
    parser.add_argument("-p", "--prefix", type = str, help = "Output Prefix for Warnings CSV")
    parser.add_argument("-f", "--format", type = str, help = "Format of Input File (overrides automatic detection)")
    parser.add_argument("--paired", action = "store_true", help = "Interleaved Paired-End Reads File")
    parser.add_argument("--strict", action = "store_true", help = "Evaluate complete file - i.e. don't output first error & exit")
    parser.add_argument("--allow_orphan_reads", action = "store_true", help = "Allow For Orphan Reads")
    argsDict = vars(parser.parse_args())

    # Set Arguments
    inputfile  = argsDict["input"]
    format = argsDict["format"]
    prefix = argsDict["prefix"]
    paired = argsDict["paired"]
    strict = argsDict["strict"]
    orphans = argsDict["allow_orphan_reads"]

    # Detect Format
    if format is None:
        format = (inputfile.split(".")[-1]).lower()

    # Get Prefix
    if prefix is None:
        prefix = inputfile.replace("." + format, "").split("/")[-1]

    # Assertions for Required Input
    assert (inputfile is not None), "No input file provided!"
    assert (format in acceptedFormats), "Invalid format: %s!" % format

    #
    # Validation
    #

    if not strict:
        warnings = {
            "unpaired_reads" : [],
            "qual_seq_len_mismatch" : [],
            "id_mismatch" : {
                "missing_1" : [],
                "missing_2" : []
            }
        }
        if paired:

            # Not Strict, Paired, No Orphans
            print("\nValidating That Sequences are Correctly Interleaved and That Sequence & Quality Scores Are of The Same Length")
            with open(inputfile, "r") as handle:
                records = SeqIO.parse(handle, format)
                resetReadframe = False
                count = 0
                for record in records:
                    seqA = record
                    seqB = next(records)
                    count += 2
                    # Are Reads Paired-End?
                    if (seqA.id[-2:] != "/1") and (seqA.id[-2:] != "/2"):
                        warnings["unpaired_reads"].append(seqA.id)
                    if (seqB.id[-2:] != "/1") and (seqB.id[-2:] != "/2"):
                        warnings["unpaired_reads"].append(seqB.id)
                    # Are Sequence Lengths & Quality Lengths the Same?
                    if (len(seqA.seq) != len(seqA.letter_annotations["phred_quality"])):
                        warnings["qual_seq_len_mismatch"].append(seqA.id)
                    if (len(seqB.seq) != len(seqB.letter_annotations["phred_quality"])):
                        warnings["qual_seq_len_mismatch"].append(seqB.id)
                    # Are Paired IDs The Same?
                    if seqA.id[0:-2] != seqB.id[0:-2]:
                        # Is Read /1 In Fact the Read /1? If Not, 1 is in Fact 2 & Missing it's /1 Pair
                        if (seqA.id[-2:] != "/1"):
                            warnings["id_mismatch"]["missing_1"].append(seqA.id)
                        # Is Read /2 In Fact the Read /2? If Not, 2 is in Fact next 1, & 1 is Missing it's /2 Pair
                        if (seqB.id[-2:] != "/2"):
                            warnings["id_mismatch"]["missing_2"].append(seqA.id)

            # Write Warnings
            with open("validation_warnings.%s.tsv" % prefix, "w") as warningfile:
                warningfile.write("id\twarning\n")
                for seqID in warnings["unpaired_reads"]:
                    warningfile.write("%s\tunpaired_reads\n" % seqID)
                for seqID in warnings["qual_seq_len_mismatch"]:
                    warningfile.write("%s\tqual_seq_len_mismatch\n" % seqID)
                for seqID in warnings["id_mismatch"]["missing_1"]:
                    warningfile.write("%s\tid_mismatch_missing_1\n" % seqID)
                for seqID in warnings["id_mismatch"]["missing_2"]:
                    warningfile.write("%s\tid_mismatch_missing_2\n" % seqID)

            # Print Results
            print("\nParsed %d Paired End Reads and Found:" % count)
            print("  %d Unpaired Read(s)" % len(warnings["unpaired_reads"]))
            print("    %d Read(s) Are Missing Their /1 Mate" % len(warnings["id_mismatch"]["missing_1"]))
            print("    %d Read(s) Are Missing Their /2 Mate" % len(warnings["id_mismatch"]["missing_2"]))
            print("  %d Read(s) With a Sequence / Quality String Length Mismatch" % len(warnings["qual_seq_len_mismatch"]))
            print("\nWarnings Written to: validation_warnings.%s.tsv" % prefix)

        else:

            # Not Strict, Single
            print("\nValidating That Single-End Sequence & Quality Scores Are of The Same Length")
            with open(inputfile, "r") as handle:
                records = SeqIO.parse(handle, format)
                count = 0
                for record in records:
                    # Are Sequence Lengths & Quality Lengths the Same?
                    assert (len(record.seq) == len(record.letter_annotations["phred_quality"])), "Error: Sequence & Quality Lengths Do Not Match: (%s)" % record.id
                    count += 1

            # Write Warnings
            with open("validation_warnings.%s.tsv" % prefix, "w") as warningfile:
                warningfile.write("id\twarning\n")
                for seqID in warnings["qual_seq_len_mismatch"]:
                    warningfile.write("%s\tqual_seq_len_mismatch\n" % seqID)
            
            # Print Results
            print("\nParsed %d Single End Reads and Found:" % count)
            print("\t%d Read(s) With a Sequence / Quality String Length Mismatch" % len(warnings["qual_seq_len_mismatch"]))
            print("\nWarnings Written to: validation_warnings.%s.tsv" % prefix)

    else:
        
        # Strict, Paired, No Orphans
        if paired:

            print("\nValidating That Sequences are Correctly Interleaved and That Sequence & Quality Scores Are of The Same Length")
            if orphans:

                with open(inputfile, "r") as handle:
                    records = SeqIO.parse(handle, format)
                    for record in records:
                        seqA = record
                        seqB = next(records)
                        # Are Reads Paired-End?
                        assert (seqA.id[-2:] == "/1") or (seqA.id[-2:] == "/2"), "Error: Sequence ID (%s) Does Not Indicate Paired Data" % seqA.id
                        assert (seqB.id[-2:] == "/1") or (seqB.id[-2:] == "/2"), "Error: Sequence ID (%s) Does Not Indicate Paired Data" % seqB.id
                        # Are Sequence Lengths & Quality Lengths the Same?
                        assert (len(seqA.seq) == len(seqA.letter_annotations["phred_quality"])), "Error: Sequence & Quality Lengths Do Not Match: (%s)" % seqA.id
                        assert (len(seqB.seq) == len(seqB.letter_annotations["phred_quality"])), "Error: Sequence & Quality Lengths Do Not Match: (%s)" % seqB.id
            
            else:

                with open(inputfile, "r") as handle:
                    records = SeqIO.parse(handle, format)
                    for record in records:
                        seqA = record
                        seqB = next(records)
                        # Are Reads Paired-End?
                        assert (seqA.id[-2:] == "/1") or (seqA.id[-2:] == "/2"), "Error: Sequence ID (%s) Does Not Indicate Paired Data" % seqA.id
                        assert (seqB.id[-2:] == "/1") or (seqB.id[-2:] == "/2"), "Error: Sequence ID (%s) Does Not Indicate Paired Data" % seqB.id
                        # Are Paired IDs The Same?
                        if seqA.id[0:-2] != seqB.id[0:-2]:
                            # Is Read /1 In Fact the Read /1? If Not, 1 is in Fact 2 & Missing it's /1 Pair
                            assert (seqA.id[-2:] == "/1"), "Orphan Read Found (Missing /1 of Pair): %s" % seqA.id
                            # Is Read /2 In Fact the Read /2? If Not, 2 is in Fact next 1, & 1 is Missing it's /2 Pair
                            assert (seqB.id[-2:] == "/2"), "Orphan Read Found (Missing /2 of Pair): %s" % seqA.id
                        # Are Sequence Lengths & Quality Lengths the Same?
                        assert (len(seqA.seq) == len(seqA.letter_annotations["phred_quality"])), "Error: Sequence & Quality Lengths Do Not Match: (%s)" % seqA.id
                        assert (len(seqB.seq) == len(seqB.letter_annotations["phred_quality"])), "Error: Sequence & Quality Lengths Do Not Match: (%s)" % seqB.id
        
        else:

            # Strict, Single
            print("\nValidating That Single-End Sequence & Quality Scores Are of The Same Length")
            with open(inputfile, "r") as handle:
                records = SeqIO.parse(handle, format)
                for record in records:
                    # Are Sequence Lengths & Quality Lengths the Same?
                    assert (len(record.seq) == len(record.letter_annotations["phred_quality"])), "Error: Sequence & Quality Lengths Do Not Match: (%s)" % record.id

    print("\nDone")

else:

    pass
