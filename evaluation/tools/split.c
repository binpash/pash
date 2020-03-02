#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1e6

void SplitInput(char* input, char* output1, char* output2, int batchSize) {
    FILE* inputFile = fopen(input, "r");
    if (!inputFile) {
        return;
    }

    FILE* outputFile1 = fopen(output1, "w");
    FILE* outputFile2 = fopen(output2, "w");
    if (!(outputFile1 && outputFile2)) {
        return;
    }
    FILE* outputFile;

    char* inputBuffer = NULL;
    unsigned int readLines = 0;

    while (getline(&inputBuffer, 0, inputFile) > 0) {
        outputFile = ++readLines >= batchSize ? outputFile2 : outputFile1;
        fputs(inputBuffer, outputFile);
    }

    fclose(inputFile);
    fclose(outputFile1);
    fclose(outputFile2);
    free(inputBuffer);
}

int main(int argc, char* argv[]) {
    // arg#1 -> input file name
    // arg#2 -> output file name 1
    // arg#3 -> output file name 2
    // arg#4 -> batch_size
    if (argc < 5) {
        printf("ERROR: missing input!");
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

    return 0;
}
