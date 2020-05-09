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

void NextFile(FILE** currentFile, char* outputFileNames[], unsigned int numOutputFiles) {
  static int current = -1;

  if (!currentFile) {
    PRINTDBG("%s: Invalid file pointer, aborting\n", __func__);
    exit(1);
  }

  if ((*currentFile)) {
    if (current >= 0) {
      PRINTDBG("%s: Will close %s output file\n", __func__, outputFileNames[current]);
    }
    fclose(*currentFile);
  }

  if (++current == numOutputFiles) {
    return;
  }

  PRINTDBG("%s: Will open %s output file\n", __func__, outputFileNames[current]);
  FILE* nextFile = fopen(outputFileNames[current], "w");
  if (!nextFile) {
    perror(LOC);
    exit(1);
  }

  PRINTDBG("%s: Successfully opened %s output file\n", __func__, outputFileNames[current]);
  *currentFile = nextFile;
}

void SplitInput(char* input, int batchSize, char* outputFileNames[], unsigned int numOutputFiles) {
  PRINTDBG("%s: will split input\n", __func__);
  FILE* outputFile = NULL;
  NextFile(&outputFile, outputFileNames, numOutputFiles);
  if (!outputFile) {
    PRINTDBG("%s: No output file in list, quitting\n", __func__);
    return;
  }

  FILE* inputFile = fopen(input, "r");
  if (!inputFile) {
    perror(LOC);
    exit(1);
  }
  PRINTDBG("%s: Opened input file %s\n", __func__, input);

  char* inputBuffer = NULL;
  unsigned int readLines = 0;

  size_t len = 0;
  while (getline(&inputBuffer, &len, inputFile) > 0) {
    if (++readLines == batchSize) {
      readLines = 0;
      NextFile(&outputFile, outputFileNames, numOutputFiles);
      if (!outputFile) {
        PRINTDBG("%s: Ran out of output files\n", __func__);
        break;
      }
      PRINTDBG("%s: successfully opened next output file\n", __func__);
    }
    fputs(inputBuffer, outputFile);
  }

  // need to exhaust our output files list
  while (outputFile) {
    NextFile(&outputFile, outputFileNames, numOutputFiles);
  }

  PRINTDBG("%s: Done splitting input %s, will clean up\n", __func__, input);
  fclose(inputFile);
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
