#include "r_split.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>
#include <sys/types.h>

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
        args[0] = "./agg_wc"; // Path to the `agg_wc` binary

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

        if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
        {
            printf("Aggregation for wc completed successfully.\n");
        }
        else
        {
            fprintf(stderr, "Aggregation failed.\n");
        }
    }
}

void inGate(char *inputFileBaseName, char *inputKeyBaseName, char *outputFileName, char *signalFileName)
{
    int version = 0;
    // open input file
    char *inputFileName = combineStrings(inputFileBaseName, "_data_v", version);
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
    int fd = open(signalFileName,  O_RDONLY | O_NONBLOCK);
    if (fd < 0)
    {
        perror("Error opening FIFO file with non-blocking mode");
        free(signalFileName);
        exit(1);
    }

    pid_t data_pid = -1;
    pid_t sig_pid = -1;

    // signal data buffer
    char signal[BUFLEN] = {0}; // Shared buffer for signal data
    ssize_t bytesRead = 0;

    size_t bufLen = BUFLEN;
    char *buffer = malloc(bufLen + 1);
    ssize_t lineLen = 0;
    int id=0;

    lineLen = getline(&buffer, &bufLen, inputFile);
    while (!feof(inputFile))
    {
        bytesRead = read(fd, signal, strlen(TIMEOUTSIG));
        if (bytesRead > 0)
        {
            signal[bytesRead] = '\0'; // Null-terminate
            printf("[ingate stateful] Signal received: %s\n", signal);

            if (strcmp(signal, TIMEOUTSIG) == 0)
            {
                printf("[ingate stateful] Timeout SIG received. Processing timeout...\n");

                /* update output file */           
                if (tempFilesCount > 0)
                {
                    fclose(outputFile);
                    printf("[FD Check] Temp file closed: %d\n", fileno(outputFile));
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
                printf("[ingate-stateful] updated outputfile fd: %d\n", fileno(outputFile));

                /* finish reading the current input file and write to the temp file */
                // lineLen = getline(&buffer, &bufLen, inputFile);
                // printf("[ingate-stateful, timed out] Reading line %ld from fd %d to fd %d\n", id, fileno(inputFile), fileno(outputFile));
                while (!feof(inputFile))
                {
                    if (fwrite(buffer, 1, lineLen, outputFile) != lineLen)
                    {
                        fprintf(stderr, "Error writing to output file: %s\n", strerror(errno));
                        free(buffer);
                        exit(EXIT_FAILURE);
                    }
                    if (outputFile && fflush(outputFile) != 0)
                    {
                        fprintf(stderr, "Error flushing file: %s (errno: %d)\n", strerror(errno), errno);
                        perror("ERror flushing \n");
                    }
                    id += 1;
                    lineLen = getline(&buffer, &bufLen, inputFile);
                }
                // printf("[ingate-stateful, timed out] Finished writing line %ld\n", id);
                safeClose(inputFile);
                inputFile = NULL;

                char *dataKey = combineStrings(inputKeyBaseName, "_v", version);
                inputFileName = combineStrings(inputFileBaseName, "_v", version);
                data_pid = createAndOpenInputFIFO(inputFileName, dataKey, &inputFile, 0666, "r");
                printf("[ingate stateful] updated input file fd: %d\n", fileno(inputFile));

                free(dataKey);

                lineLen = getline(&buffer, &bufLen, inputFile);
                if (feof(inputFile)) {
                    printf("[ingate stateful] eof break\n");
                    break;
                }
            } // if timeout 
        } // if byte > 0

        // Normal operation: write data to output file
        // if (tempFilesCount > 0)
        //     printf("[ingate stateful] Reading line %ld from fd %d to fd %d\n", id, fileno(inputFile), fileno(outputFile));

        if (fwrite(buffer, 1, lineLen, outputFile) != lineLen)
        {
            fprintf(stderr, "[ingate stateful] Error writing to %s: %s\n", tempFiles[tempFilesCount - 1], strerror(errno));
            free(buffer);
            exit(EXIT_FAILURE);
        }
        if (outputFile && fflush(outputFile) != 0)
        {
            fprintf(stderr, "Error flushing file: %s (errno: %d)\n", strerror(errno), errno);
            perror("ERror flushing \n");
        }
        // printf("[ingate stateful] Finished writing line %ld\n", id);
        id += 1;

        lineLen = getline(&buffer, &bufLen, inputFile);
    }

    int isWC = 1;
    if (isWC == 1){
        // wc
        aggregateWC(tempFiles, tempFilesCount, outputFileSaved);
    } else {
        // sort
        aggregateFiles(tempFiles, tempFilesCount, outputFileSaved); 
    }
    printf("[ingate-stateful] Finish aggregating files\n");

    for (int i = 0; i < tempFilesCount; i++) {
        free(tempFiles[i]);
    }

    // Free the inputFileName and signalFileName if allocated dynamically
    free(inputFileName);
    free(signalFileName);
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
    char *signalFileName = argv[4];

    inGate(inputFileBaseName, inputKeyBaseName, outputFileName, signalFileName);

    // clean up
    return 0;
}