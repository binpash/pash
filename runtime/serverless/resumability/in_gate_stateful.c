#include "r_split.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>
#include <sys/types.h>

#define BUFFER_SIZE 100000

void aggregateFiles(char **temp_files, int numFiles, FILE *outputFile)
{
    if (numFiles == 0)
    {
        printf("[aggregateFiles] No timeout\n");
        return;
    } else if (temp_files == NULL || outputFile == NULL)
    {
        fprintf(stderr, "[aggregateFiles] Invalid arguments\n");
        return;
    }

    pid_t pid = fork();

    if (pid < 0)
    {
        perror("Failed to fork");
        exit(EXIT_FAILURE);
    }
    else if (pid == 0)
    {                                           
        if (dup2(fileno(outputFile), STDOUT_FILENO) < 0) {
            perror("[aggregateFiles] dup2 failed");
            exit(EXIT_FAILURE);
        }
        fclose(outputFile); // Close the FILE* as it's no longer needed

        // Prepare arguments for execvp
        char **args = malloc((numFiles + 3) * sizeof(char *));
        if (args == NULL)
        {
            perror("Failed to allocate memory for arguments");
            exit(EXIT_FAILURE);
        }

        args[0] = "sort";
        args[1] = "-m";
        for (int i = 0; i < numFiles; i++)
        {
            args[i + 2] = temp_files[i];
        }
        args[numFiles + 2] = NULL; // Null-terminate the arguments

        // Execute the sort command
        execvp("sort", args);

        // If execvp fails
        perror("[aggregateFiles] execvp failed");
        free(args);
        exit(EXIT_FAILURE);
    }
    else
    {
        int status;
        if (waitpid(pid, &status, 0) < 0)
        {
            perror("Failed to wait for child process");
            exit(EXIT_FAILURE);
        }

        if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
        {
            printf("Aggregation completed successfully.\n");
        }
        else
        {
            fprintf(stderr, "Aggregation failed.\n");
        }
    }
}

void aggregateWC(char **temp_files, int numFiles, FILE *outputFile)
{
    // printf("[aggregateWC] Start aggregating files\n");
    if (numFiles == 0)
    {
        printf("[aggregateFiles] No files to process\n");
        return;
    }
    else if (temp_files == NULL || outputFile == NULL)
    {
        fprintf(stderr, "[aggregateFiles] Invalid arguments\n");
        return;
    }

    pid_t pid = fork();

    if (pid < 0)
    {
        perror("Failed to fork");
        exit(EXIT_FAILURE);
    }
    else if (pid == 0)
    {
        // Child process: Redirect stdout to the outputFile
        int fd = fileno(outputFile);
        if (fd == -1)
        {
            perror("Failed to get file descriptor for output file");
            exit(EXIT_FAILURE);
        }

        // Duplicate the output file descriptor to stdout
        if (dup2(fd, STDOUT_FILENO) == -1)
        {
            perror("Failed to redirect stdout to output file");
            exit(EXIT_FAILURE);
        }

        // Prepare arguments for the `agg_wc` binary
        char *args[numFiles + 2];
        args[0] = "./runtime/agg-wc-l.sh"; // Path to the `agg_wc` binary

        for (int i = 0; i < numFiles; i++)
        {
            args[i + 1] = temp_files[i];
        }
        args[numFiles + 1] = NULL; // NULL-terminated array

        // Execute the `agg_wc` binary
        if (execvp(args[0], args) == -1)
        {
            perror("Failed to execute agg_wc");
            exit(EXIT_FAILURE);
        }

        exit(EXIT_SUCCESS);
    }
    else
    {
        // Parent process: Wait for the child to complete
        int status;
        if (waitpid(pid, &status, 0) < 0)
        {
            perror("Failed to wait for child process");
            exit(EXIT_FAILURE);
        }

        if (!(WIFEXITED(status) && WEXITSTATUS(status) == 0)) {
            fprintf(stderr, "Aggregation failed.\n");
        }
    }
}



void inGate(char *inputFileBaseName, char *inputKeyBaseName, char *outputFileName, char *lambdaScriptID)
{
    int version = 0;

    // open input file
    char *inputFileName = combineStrings(inputFileBaseName, "_v", version);
    FILE *inputFile = fopen(inputFileName, "r");
    if (!inputFile)
    {
        perror("Error opening input file");
        exit(1);
    }

    // open output file
    FILE *outputFile = fopen(outputFileName, "w");
    FILE *outputFileSaved = outputFile;
    if (!outputFile)
    {
        perror("Error opening output file");
        exit(1);
    }
    // assume filename length is less than 100
    int tempFilesCount = 0;
    char *tempFiles[100];

    // open signal file: non blocking mode
    char *signalFileName = combineStrings(lambdaScriptID, "_sig_", 0);
    int fd = open(signalFileName,  O_RDONLY | O_NONBLOCK);
    while (fd < 0) // ideally check the type of error and only retry on ENOENT
    {
        // perror("Error opening FIFO file with non-blocking mode");
        // free(signalFileName);
        // exit(1);
        fd = open(signalFileName,  O_RDONLY | O_NONBLOCK);
    }

    pid_t data_pid = -1;

    // signal data buffer
    char signal[BUFLEN] = {0}; // Shared buffer for signal data
    ssize_t bytesRead = 0;

    size_t bufLen = BUFFER_SIZE;
    char *buffer = malloc(bufLen + 1);
    ssize_t len = 0;

    // printf("[ingate stateful %s] Start reading input file\n", inputKeyBaseName);
    // log_with_timestamp("[ingate stateful %s] Start reading input file", inputKeyBaseName);

    while ((len = fread(buffer, 1, bufLen, inputFile)) > 0 )
    {
        // printf("[ingate stateful %s] Read %ld bytes\n", inputKeyBaseName, len);
        // check signal
        bytesRead = read(fd, signal, strlen(TIMEOUTSIG));
        if (bytesRead > 0)
        {
            signal[bytesRead] = '\0'; // Null-terminate

            if (strcmp(signal, TIMEOUTSIG) == 0)
            {
                // first time recv timeout, redirect the output then write
                if (tempFilesCount == 0) 
                {
                    tempFilesCount += 1;
                    char *tempFileName = combineStrings(outputFileName, "_temp", tempFilesCount);
                    tempFiles[tempFilesCount - 1] = tempFileName;
                    outputFile = fopen(tempFileName, "w+"); // exclusive access
                    if (!outputFile) {
                        perror("Failed to open temporary file");
                        exit(EXIT_FAILURE);
                    }
                    /* finish reading the current input file and write to the temp file */
                    safeWriteWithFlush(buffer, 1, len, outputFile);
                    while (!feof(inputFile))
                    {
                        len = fread(buffer, 1, bufLen, inputFile);
                        safeWriteWithFlush(buffer, 1, len, outputFile);
                    }
                    // printf("[ingate-stateful, timed out] Finished writing line %ld\n", id);
                    safeClose(inputFile);
                    inputFile = NULL;
                }        
                else 
                {   
                    // not first time, read everything, then redirect to the next temp file
                    safeWriteWithFlush(buffer, 1, len, outputFile);
                    while (!feof(inputFile))
                    {
                        len = fread(buffer, 1, bufLen, inputFile);
                        safeWriteWithFlush(buffer, 1, len, outputFile);
                    }
                    // printf("[ingate-stateful, timed out] Finished writing line %ld\n", id);
                    safeClose(inputFile);
                    inputFile = NULL;

                    // close the temp file
                    fclose(outputFile);
                    // printf("[FD Check] Temp file closed: %d\n", fileno(outputFile));
                    outputFile = NULL;
                }
                tempFilesCount += 1;
                char *tempFileName = combineStrings(outputFileName, "_temp", tempFilesCount);
                tempFiles[tempFilesCount - 1] = tempFileName;
                outputFile = fopen(tempFileName, "w+"); // exclusive access
                if (!outputFile) {
                    perror("Failed to open temporary file");
                    exit(EXIT_FAILURE);
                }
                
                // modify createAndOpenInputFIFO to returna data_pid too 
                if (data_pid != -1)
                {
                    waitpid(data_pid, NULL, 0);
                }

                /* spawning new pashlib process */
                version += 1;
                char *dataKey = combineStrings(inputKeyBaseName, "_v", version);
                inputFileName = combineStrings(inputFileBaseName, "_v", version);
                data_pid = createAndOpenInputFIFO(inputFileName, dataKey, &inputFile, 0666, "r");
                // printf("[ingate stateful] updated input file fd: %d\n", fileno(inputFile));
                free(dataKey);
            }  else {
                safeWriteWithFlush(buffer, 1, len, outputFile);
            }
        } else {
            safeWriteWithFlush(buffer, 1, len, outputFile);
        }
    }

    safeClose(inputFile);

    if (tempFilesCount != 0)
    {
        int isWC = 1;
        if (isWC == 1){
            // wc
            aggregateWC(tempFiles, tempFilesCount, outputFileSaved);
        } else {
            // sort
            aggregateFiles(tempFiles, tempFilesCount, outputFileSaved); 
        }

        for (int i = 0; i < tempFilesCount; i++) {
        free(tempFiles[i]);
    }
    }
    // clean up
    
    safeClose(outputFileSaved);

    // Free the inputFileName and signalFileName if allocated dynamically
    free(inputFileName);
}

// arg#1: data input file basename
// arg#2: data input key basename
// arg#3: data output file
// arg#4: signal file
int main(int argc, char *argv[])
{
    char *inputFileBaseName = argv[1];
    char *inputKeyBaseName = argv[2];
    char *outputFileName = argv[3];
    char *lambdaSriptID = argv[4];

    inGate(inputFileBaseName, inputKeyBaseName, outputFileName, lambdaSriptID);

    // clean up
    return 0;
}