#include "r_split.h"

typedef struct block
{
  int64_t id;
  size_t blockSize;
  FILE *inputFile;
} block_t;

//prints to stdout, could be modified to write to a file
void MergeInput(char *inputFileNames[], unsigned int numInputFiles)
{
  //ring buffer, this implemntation assumes that each round over the file will have the next complete set of blocks
  //(ie in 2 blocks, block 1 will have ids 0..2..4 and block 2 will have 1..3..5)
  block_t *blockBuf = malloc(sizeof(block_t) * numInputFiles);
  size_t bufLen = BUFLEN; //buffer length, would be resized as needed
  char *buffer = malloc(bufLen + 1);
  int64_t id;
  size_t blockSize;
  bool isLast;

  for (int i = 0; i < numInputFiles; i++)
  {
    FILE *inputFile = fopen(inputFileNames[i], "r");
    if (!inputFile)
    {
      perror(LOC);
      exit(1);
    }
    blockBuf[i].inputFile = inputFile;
  }
  int bufferIdx = 0;
  int64_t nextID = 0;
  for(;;) {
    //read header
    readHeader(blockBuf[bufferIdx].inputFile, &id, &blockSize, &isLast);
    if (feof(blockBuf[bufferIdx].inputFile)) {
      PRINTDBG("r_merge: End of file");
      break;
    }
    blockBuf[bufferIdx].id = id;
    blockBuf[bufferIdx].blockSize = blockSize;
    
    size_t tot_read = 0, readSize = 0;
    // fprintf(stderr, "start reading block %ld\n", id);
    while (tot_read < blockSize) {
        readSize = MIN(bufLen, blockSize-tot_read);

        readSize = handle_reading(buffer, readSize, blockBuf[bufferIdx].inputFile);
        
        tot_read += readSize;
        safeWrite(buffer, 1, readSize, stdout);
    }
    fflush(stdout); //only need to flush after block ended
    assert(tot_read == blockSize);
    
    if (isLast) {
      bufferIdx = (bufferIdx + 1) % numInputFiles;
      nextID += 1;
    }
  };

  //assert reason while loop terminated is eof
  assert(feof(blockBuf[bufferIdx].inputFile));

  //clean up
  for (int i = 0; i < numInputFiles; i++)
  {
    fclose(blockBuf[i].inputFile);
  }
  free(blockBuf);
  free(buffer);
}

int main(int argc, char *argv[])
{
  // args#2... -> input file names

  if (argc < 2)
  {
    // TODO print usage string
    fprintf(stderr, "missing input!\n");
    exit(1);
  }

  char **inputFileNames = (char **)malloc((argc - 1) * sizeof(char *));
  for (int i = 1; i < argc; i++)
  {
    inputFileNames[i - 1] = calloc(strlen(argv[i]) + 1, sizeof(char));
    strcpy(inputFileNames[i - 1], argv[i]);
  }

  MergeInput(inputFileNames, argc - 1);

  return 0;
}
