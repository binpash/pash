#include "r_split.h"

void SplitByBytes(FILE *inputFile, int batchSize, FILE *outputFiles[], unsigned int numOutputFiles)
{
  int current = 0;
  int64_t id = 0;
  size_t len = 0;
  FILE *outputFile = outputFiles[current];

  char *buffer = malloc(batchSize + 1);

  // Do round robin copying of the input file to the output files
  // Each block has a header of "ID blockSize\n"
  while ((len = fread(buffer, 1, batchSize, inputFile)) > 0)
  {
    //write header
    writeHeader(outputFile, id, len);

    //write blocks
    fwrite(buffer, 1, len, outputFile);

    current = (current + 1) % numOutputFiles;
    outputFile = outputFiles[current];
    id += 1;
  }

  if (len < 0)
  {
    perror(LOC);
    exit(1);
  }

  //clean up
  free(buffer);
}

void SplitByLines(FILE *inputFile, int batchSize, FILE *outputFiles[], unsigned int numOutputFiles, int8_t add_header)
{
  int current = 0;
  int64_t id = 0;
  size_t len = 0, headSize = 0, restSize = 0, prevRestSize = 0, blockSize = 0, bufLen = 0;
  FILE *outputFile = outputFiles[current];

  char *buffer = malloc(batchSize + 1);
  char *incompleteLine = malloc(batchSize + 1);
  char *newLineBuffer = NULL; //only used when a newline character is not found in chunk

  // Do round robin copying of the input file to the output files
  // Each block has a header of "ID blockSize\n"
  while ((len = fread(buffer, 1, batchSize, inputFile)) > 0)
  {
    //find pivot point for head and rest
    for (int i = len - 1; i >= 0; i--)
    {
      if (buffer[i] == '\n')
      {
        headSize = i + 1;
        restSize = len - headSize;
        break;
      }
    }

    //no new line character
    if (headSize == 0)
    {
      headSize = len;
      if ((len = getline(&newLineBuffer, &bufLen, inputFile)) < 0)
      {
        //edge case to fix: can't be called if file ended
        fprintf(stderr, "r_split: getline failed");
        exit(1);
      }
      blockSize = prevRestSize + headSize + len;
      if (add_header)
        writeHeader(outputFile, id, blockSize);
      //write blocks
      if (prevRestSize)
        safeWrite(incompleteLine, 1, prevRestSize, outputFile);
      safeWrite(buffer, 1, headSize, outputFile);
      safeWriteWithFlush(newLineBuffer, 1, len, outputFile);
    }
    else
    {
      blockSize = prevRestSize + headSize;
      //write header
      if (add_header)
        writeHeader(outputFile, id, blockSize);
      //write blocks
      if (prevRestSize)
        safeWrite(incompleteLine, 1, prevRestSize, outputFile);
      safeWriteWithFlush(buffer, 1, headSize, outputFile);
      //update incompleteLine to the current block
      memcpy(incompleteLine, buffer + headSize, restSize);
    }
    // fflush(outputFile);
    current = (current + 1) % numOutputFiles;
    outputFile = outputFiles[current];
    prevRestSize = restSize;
    headSize = restSize = 0;
    id += 1;
  }

  if (len < 0)
  {
    perror(LOC);
    exit(1);
  }

  if (prevRestSize > 0)
  {
    if (add_header)
      writeHeader(outputFile, id, prevRestSize);
    safeWriteWithFlush(incompleteLine, 1, prevRestSize, outputFile);
  }

  //clean up
  free(buffer);
  free(incompleteLine);
  if (newLineBuffer)
    free(newLineBuffer);
}

void SplitByLinesRaw(FILE *inputFile, int batchSize, FILE *outputFiles[], unsigned int numOutputFiles)
{
  int current = 0, len = 0;
  size_t bufLen = 0;
  FILE* outputFile = outputFiles[current];

  char* buffer = NULL;

  // Do round robin copying of the input file to the output files without any headers
  while ((len = getline(&buffer, &bufLen, inputFile)) > 0) {
      safeWrite(buffer, 1, len, outputFile);
      current = (current + 1) % numOutputFiles;
      outputFile = outputFiles[current];
  }
  //clean up
  free(buffer);

}



void SplitInput(char *input, int batchSize, char *outputFileNames[], unsigned int numOutputFiles, int8_t useBytes, int8_t raw)
{
  PRINTDBG("%s: will split input\n", __func__);
  //TODO: find better way?
  FILE **outputFiles = malloc(sizeof(FILE *) * numOutputFiles);
  for (int i = 0; i < numOutputFiles; i++)
  {
    outputFiles[i] = fopen(outputFileNames[i], "w");
    if (!outputFiles[i])
    {
      perror(LOC);
      exit(1);
    }
  }

  FILE *inputFile = fopen(input, "r");
  if (!inputFile)
  {
    perror(LOC);
    exit(1);
  }
  PRINTDBG("%s: Opened input file %s\n", __func__, input);

  if (raw)
  {
    // SplitByLinesRaw(inputFile, batchSize, outputFiles, numOutputFiles);
    SplitByLines(inputFile, batchSize, outputFiles, numOutputFiles, 0);
  }
  else if (useBytes)
  {
    //if batchSize isn't set we can approximate it
    if (batchSize == 0)
    {
      int inputfd = fileno(inputFile);
      struct stat buf;
      fstat(inputfd, &buf);
      size_t inputSize = buf.st_size;
      batchSize = MIN(inputSize / MINCHUNKS, CHUNKSIZE); //autotune this better?
    }
    SplitByBytes(inputFile, batchSize, outputFiles, numOutputFiles);
  }
  else
  {
    SplitByLines(inputFile, batchSize, outputFiles, numOutputFiles, 1);
  }

  PRINTDBG("%s: Done splitting input %s, will clean up\n", __func__, input);
  fclose(inputFile);
  // need to close all output files
  for (int i = 0; i < numOutputFiles; i++)
  {
    fclose(outputFiles[i]);
  }
  free(outputFiles);
}

int main(int argc, char *argv[])
{
  // arg#1 -> input file name
  // arg#2 -> batch_size
  // args#3... -> output file names
  // flags: -b to use bytes (batch_size will be exact number of bytes instead of approximating to the closest line)
  if (argc < 4)
  {
    // TODO print usage string
    fprintf(stderr, "missing input!\n");
    exit(1);
  }
  int8_t useBytes = 0, offset = 0, raw = 0;
  size_t batchSize = 0;
  char **outputFileNames = NULL;
  char *inputFileName = NULL;
  for (int i = 1; i < argc; i++)
  {
    //check flags
    if (strcmp(argv[i], "-b") == 0)
    {
      useBytes = 1;
      offset += 1;
      continue;
    }

    if (strcmp(argv[i], "-r") == 0)
    {
      raw = 1;
      offset += 1;
      continue;
    }

    int truePos = i - offset; //tracks true position without any added flags
    if (truePos == 1)
    {
      inputFileName = calloc(strlen(argv[i]) + 1, sizeof(char));
      strcpy(inputFileName, argv[i]);
    }

    if (truePos == 2)
      batchSize = atoi(argv[i]); //if 0 r_split default to approximate value

    if (truePos == 3)
    {
      outputFileNames = argv + i;
    }
  }

  SplitInput(inputFileName, batchSize, outputFileNames, argc - offset - 3, useBytes, raw);

  PRINTDBG("SplitInput is done\n");
  return 0;
}
