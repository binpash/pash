#! /home/users/cordier/.linuxbrew/bin/python3

if __name__ == "__main__":

    # Imports
    import argparse
    import random
    import os
    # Library Imports
    from Bio import SeqIO

    # 
    # Parse Arguments
    # 

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "Input Paired-End FastQ File")
    parser.add_argument("-o", "--outputprefix", type = str, help = "Prefix for Output FastQ File (i.e. prefix.shuffled.fastq)")
    parser.add_argument("-s", "--chunksize", type = int, help = "Size of Chunk (i.e. n Pairs of Reads) to Shuffle at a Time. Default = 50000")
    argsDict = vars(parser.parse_args())

    fastq = argsDict["input"]
    prefix = argsDict["outputprefix"]
    chunksize = argsDict["chunksize"]

    # Assertions for Required Input
    assert (fastq is not None), "No FastQ input provided!"

    #  If No Prefix, Use Same as FASTQ
    if prefix is None: 
        prefix = ".".join(fastq.split(".")[0:-1])

    # If No Chunk Size, Set Default
    if chunksize is None:
        chunksize = 100000
    assert (chunksize % 2 == 0), "Chunksize must be even!"
    
    # 
    # Conversion
    #

    # Open FastQ, Chunk File, Shuffle, & Write Chunk to Output
    with open(fastq, "r") as handle:
        records = SeqIO.parse(handle, "fastq")
        chunk = []
        chunks, nrecords = 0, 0
        for record in records:
            try:
                chunk.append([record, next(records)])
                if len(chunk) > chunksize:
                    with open(prefix + ".shuffled.chunk_%d.fastq" % chunks, "w") as output:
                        random.shuffle(chunk) # Inplace
                        flattened = [item for sublist in chunk for item in sublist]
                        SeqIO.write(flattened, output, "fastq")
                        nrecords += len(flattened)
                        print("Shuffle FastQ: Chunk Written to %s.shuffled.chunk_%d.fastq (%d Records)" % (prefix, chunks, nrecords))
                        chunks += 1
                        chunk = []
            except:
                with open(prefix + ".shuffled.chunk_%d.fastq" % chunks, "w") as output:
                    random.shuffle(chunk) # Inplace
                    flattened = [item for sublist in chunk for item in sublist]
                    SeqIO.write(flattened, output, "fastq")
                    nrecords += len(flattened)
                    print("Shuffle FastQ: Chunk Written to %s.shuffled.chunk_%d.fastq (%d Records)" % (prefix, chunks, nrecords))
                    chunks += 1
                    chunk = []
        # Write Tail
        try:
            with open(prefix + ".shuffled.chunk_%d.fastq" % chunks, "w") as output:
                random.shuffle(chunk) # Inplace
                flattened = [item for sublist in chunk for item in sublist]
                SeqIO.write(flattened, output, "fastq")
                nrecords += len(flattened)
                print("Shuffle FastQ: Chunk Written to %s.shuffled.chunk_%d.fastq (%d Records)" % (prefix, chunks, nrecords))
                chunks += 1
        except:
            pass

    print("Shuffle FastQ: Shuffling Chunks & Reassembling FastQ")
    # Write Chunk Files in Random Order
    with open(prefix + ".shuffled.fastq", "w") as output, open(prefix + ".shuffle_order.txt", "w") as order:
        i = 0
        for chunk in random.sample(range(chunks), chunks):
            order.write(chunk)
            with open(prefix + ".shuffled.chunk_%d.fastq" % chunk, "r") as chunkfile:
                records = SeqIO.parse(chunkfile, "fastq")
                SeqIO.write(records, output, "fastq")
            # Delete Chunk
            os.unlink(prefix + ".shuffled.chunk_%d.fastq" % chunk)
            i += 1
            print("Shuffle FastQ: Chunk %d Written (%d / %d)" % (chunk, i, chunks))



else:

    pass
