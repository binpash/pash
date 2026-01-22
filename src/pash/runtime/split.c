#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1e6

#define LOC __FILE__ // TODO expand

#ifdef DEBUG
#define PRINTDBG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define PRINTDBG(fmt, ...)
#endif

void SplitInput(char* input, int batchSize, char* outputFileNames[], unsigned int numOutputFiles) {
  PRINTDBG("%s: will split input\n", __func__);

  PRINTDBG("%s: Openning output files\n", __func__);
  int current = 0;
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
  FILE* outputFile = outputFiles[current];
  
  FILE *inputFile = fopen(input, "r");
  if (!inputFile)
  {
    perror(LOC);
    exit(1);
  }
  PRINTDBG("%s: Opened input file %s\n", __func__, input);

  char* inputBuffer = NULL;
  unsigned int readLines = 0;

  size_t len = 0;
  while (getline(&inputBuffer, &len, inputFile) > 0 && !ferror(outputFile)) {
    if (++readLines == batchSize && current < numOutputFiles - 1) {
      readLines = 0;
      fclose(outputFile);
      current += 1;
      outputFile = outputFiles[current];
    }
    fputs(inputBuffer, outputFile);
  }

  PRINTDBG("%s: Done splitting input %s, will clean up\n", __func__, input);
  fclose(inputFile);
  // need to close the last output file
  fclose(outputFile);
  free(inputBuffer);
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

  int batchSize = atoi(argv[2]);

  char** outputFileNames = (char **) malloc((argc - 3) * sizeof(char *));
  for (int i = 3; i < argc; i++) {
    outputFileNames[i - 3] = calloc(strlen(argv[i]) + 1, sizeof(char));
    strcpy(outputFileNames[i - 3], argv[i]);
  }

  SplitInput(inputFileName, batchSize, outputFileNames, argc - 3);

  PRINTDBG("SplitInput is done\n");
  return 0;
}
