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

void SplitInput(char* input, char* output1, char* output2, int batchSize) {
  PRINTDBG("will split input\n");
  FILE* inputFile = fopen(input, "r");
  if (!inputFile) {
    perror(LOC);
    exit(1);
  }
  PRINTDBG("opened input file %s\n", input);

  PRINTDBG("will open outputFile1 from %s \n", output1);
  FILE* outputFile1 = fopen(output1, "w");
  if (!outputFile1) {
    perror(LOC);
    exit(1);
  }
  PRINTDBG("opened output file %s\n", output1);
  FILE* outputFile = outputFile1;

  char* inputBuffer = NULL;
  unsigned int readLines = 0;

  size_t len = 0;
  while (getline(&inputBuffer, &len, inputFile) > 0) {
    if (++readLines == batchSize) {
      fclose(outputFile);
      if (!(outputFile = fopen(output2, "w"))) {
        perror(LOC);
        exit(1);
      }
      PRINTDBG("successfully opened second output file %s\n", output2);
    }
    fputs(inputBuffer, outputFile);
  }


  PRINTDBG("done reading, will close all open files from %s, %s, %s\n", input, output1, output2);
  fclose(inputFile);
  fclose(outputFile);
  free(inputBuffer);
}

int main(int argc, char* argv[]) {
  // arg#1 -> input file name
  // arg#2 -> output file name 1
  // arg#3 -> output file name 2
  // arg#4 -> batch_size
  if (argc < 5) {
    fprintf(stderr, "missing input!\n");
    exit(1);
  }

  char* inputFileName = calloc(strlen(argv[1]) + 1, sizeof(char));
  strcpy(inputFileName, argv[1]);

  char* outputFileName1 = calloc(strlen(argv[2]) + 1, sizeof(char));
  strcpy(outputFileName1, argv[2]);
  char* outputFileName2 = calloc(strlen(argv[3]) + 1, sizeof(char));
  strcpy(outputFileName2, argv[3]);

  int batchSize = atoi(argv[4]);

  SplitInput(inputFileName, outputFileName1, outputFileName2, batchSize);

  PRINTDBG("SplitInput is done, whatever that's worth\n");
  return 0;
}
