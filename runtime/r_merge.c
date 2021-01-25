#include "r_split.h"


typedef struct block {
  int64_t id;
  size_t blockSize;
  FILE* inputFile;
} block_t;

//prints to stdout, could be modified to write to a file
void MergeInput(char* inputFileNames[], unsigned int numInputFiles) {
    //ring buffer, this implemntation assumes that each round over the file will have the next complete set of blocks 
    //(ie in 2 blocks, block 1 will have ids 0..2..4 and block 2 will have 1..3..5) 
    block_t* blockBuf = malloc(sizeof(block_t) * numInputFiles);
    size_t bufLen = CHUNKSIZE; //buffer length, would be resized as needed
    char* buffer = malloc(bufLen + 1);
    int64_t id;
    size_t blockSize;
    
    for(int i = 0; i < numInputFiles; i++) {
        FILE* inputFile = fopen(inputFileNames[i], "r");
        if (!inputFile) {
            perror(LOC);
            exit(1);
        }
        readHeader(inputFile, &id, &blockSize);
        //ordering the buffer based on id (could be unnecessary)
        blockBuf[id].inputFile = inputFile;
        blockBuf[id].id = id;
        blockBuf[id].blockSize = blockSize;
    }

    int bufferIdx = 0;
    int64_t nextID = 0;
    while (blockBuf[bufferIdx].id == nextID) {
      if(blockBuf[bufferIdx].blockSize > bufLen) {
        bufLen = blockBuf[bufferIdx].blockSize;
        buffer = realloc(buffer, bufLen + 1);
      }
      if (fread(buffer, 1, blockBuf[bufferIdx].blockSize, blockBuf[bufferIdx].inputFile) != blockBuf[bufferIdx].blockSize) {
        fprintf(stderr, "There is a problem with the merge files head data");
        exit(1);
      }
      
      fwrite(buffer, 1, blockBuf[bufferIdx].blockSize, stdout);
      
      if (feof(blockBuf[bufferIdx].inputFile)) {
        blockBuf[bufferIdx].id = -1;
      }
      
      //update contents to reflect the next block
      readHeader(blockBuf[bufferIdx].inputFile, &id, &blockSize);
     
      blockBuf[bufferIdx].id = id;
      blockBuf[bufferIdx].blockSize = blockSize;

      bufferIdx = (bufferIdx + 1) % numInputFiles;
      nextID += 1;
    }

    //clean up
    for(int i = 0; i < numInputFiles; i++) {
      fclose(blockBuf[i].inputFile);
    }
    free(blockBuf);
    free(buffer);
}

int main(int argc, char* argv[]) {
  // args#2... -> input file names
  
  if (argc < 2) {
    // TODO print usage string
    fprintf(stderr, "missing input!\n");
    exit(1);
  }
  
  char** inputFileNames = (char **) malloc((argc - 1) * sizeof(char *));
  for (int i = 1; i < argc; i++) {
    inputFileNames[i - 1] = calloc(strlen(argv[i]) + 1, sizeof(char));
    strcpy(inputFileNames[i - 1], argv[i]);
  }
  
  MergeInput(inputFileNames, argc - 1);
  
  return 0;
}
