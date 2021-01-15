#include "r_split.h"

#ifdef DEBUG
#define PRINTDBG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define PRINTDBG(fmt, ...)
#endif

//TODO: batchs should always end with a new line(if batchSize is too small now, the invariant is broken)
void SplitInput(char* input, int batchSize, char* outputFileNames[], unsigned int numOutputFiles) {
  PRINTDBG("%s: will split input\n", __func__);
  //TODO: find better way?
  FILE* outputFiles[NUMOUTFILESLIMIT] = {NULL};
  for(int i = 0; i < numOutputFiles; i++) {
    outputFiles[i] = fopen(outputFileNames[i], "w");
    if (!outputFiles[i]) {
      perror(LOC);
      exit(1);
    }
  }
  int current = 0;

  FILE* inputFile = fopen(input, "r");
  if (!inputFile) {
    perror(LOC);
    exit(1);
  }
  PRINTDBG("%s: Opened input file %s\n", __func__, input);
  
  size_t len = 0, totalbytes = 0, headSize = 0, restSize = 0, prevRestSize = 0, writeSize = 0;
  FILE* outputFile = outputFiles[current];
  size_t readSize = batchSize;

  //if size isn't set we can approximate it
  if (!readSize) {
    int inputfd = fileno(inputFile);
    struct stat buf;
    fstat(inputfd, &buf);
    size_t inputSize = buf.st_size;
    readSize = MIN(inputSize/MINCHUNKS, CHUNKSIZE); //autotune this better or replace with batchSize?
  }

  //allocate both buffers at the same time
  char* buffer = malloc(2*readSize + 2);
  char* incompleteLine = buffer + readSize + 1;
  int blockID = 0;

  

  // Research: find if there is benefit to fflush
  // Do round robin copying of the input file to the output files
  // Each block has a header of one number representing the size of the next chunk block
  while ((len = fread(buffer, 1, readSize, inputFile)) > 0) {

    //find pivot point for head and rest
    for (int i = len - 1; i >= 0; --i) {
      if (buffer[i] == '\n') {
        headSize = i + 1;
        restSize = len - headSize;
        break;
      }
    }

    if (headSize == 0) {
      headSize = len;
    }
    assert(len == (headSize + restSize));

    //write header
    writeSize = prevRestSize + headSize;
    fprintf(outputFile, "%d %lu\n", blockID++, writeSize); //try writing the actual bytes instead (true stream)
    
    //write blocks
    if (prevRestSize)
      fwrite(incompleteLine, 1, prevRestSize, outputFile);
    fwrite(buffer, 1, headSize, outputFile);
    
    //update incompleteLine to the current block
    memcpy(incompleteLine, buffer + headSize , restSize);

    current = (current + 1) % numOutputFiles;
    outputFile = outputFiles[current];
    prevRestSize = restSize;
    headSize = restSize = 0;

    #ifdef DEBUG
      totalbytes += writeSize;
    #endif
  }

  if (readSize < 0) {
    perror(LOC);
    exit(1);
  }

  PRINTDBG("%s: Done splitting input %s, will clean up\n", __func__, input);
  fclose(inputFile);
  // need to close all output files
  for(int i = 0; i < numOutputFiles; i++) {
    fclose(outputFiles[i]);
  };
  free(buffer);
}

int main(int argc, char* argv[]) {
  // arg#1 -> input file name
  // arg#2 -> batch_size
  // args#3... -> output file names
  if (argc < 4) {
    // TODO print usage string
    fprintf(stderr, "missing input!\n");
    exit(1);
  }

  char* inputFileName = calloc(strlen(argv[1]) + 1, sizeof(char));
  strcpy(inputFileName, argv[1]);

  int batchSize = atoi(argv[2]); //if 0 r_split default to approximate value

  char** outputFileNames = (char **) malloc((argc - 3) * sizeof(char *));
  for (int i = 3; i < argc; i++) {
    outputFileNames[i - 3] = calloc(strlen(argv[i]) + 1, sizeof(char));
    strcpy(outputFileNames[i - 3], argv[i]);
  }

  SplitInput(inputFileName, batchSize, outputFileNames, argc - 3);

  PRINTDBG("SplitInput is done\n");
  return 0;
}
