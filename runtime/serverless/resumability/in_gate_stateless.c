#include "r_split.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>
#include <sys/types.h>

void inGate(char *inputFileBaseName, char *inputKeyBaseName, char *outputFileName, char *lambda_script_id)
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
    if (!outputFile)
    {
        perror("Error opening output file");
        exit(1);
    }

    // open signal file: non blocking mode
    char *signalFileName = combineStrings(lambda_script_id, "_sig_", 0);
    // printf("[ingate stateless] Opening signal file: %s\n", signalFileName);
    // if (access(signalFileName, F_OK) == -1)
    // {
    //     printf("[ingate stateless] File does not exist. Creating %s.\n", filename);
    //     if (mkfifo(filename, fifo_mode) < 0)
    //     {
    //     perror("mkfifo");
    //     free(filename);
    //     exit(1);
    //     }
    // }
    int fd = open(signalFileName,  O_RDONLY | O_NONBLOCK);
    while (fd < 0)
    {
        // perror("Error opening FIFO file with non-blocking mode");
        // free(signalFileName);
        // exit(1);
        fd = open(signalFileName,  O_RDONLY | O_NONBLOCK);
    }
    // printf("[ingate stateless] Opened signal file: %s\n", signalFileName);

    pid_t data_pid = -1;

    // signal data buffer
    char signal[BUFLEN] = {0}; // Shared buffer for signal data
    ssize_t bytesRead = 0;

    size_t bufLen = BUFLEN;
    char *buffer = malloc(bufLen + 1);
    int64_t id;
    size_t blockSize;
    bool isLast;

    printf("[ingate stateless %s] Start reading input file\n", inputKeyBaseName);

    readHeader(inputFile, &id, &blockSize, &isLast);
    while (!feof(inputFile))
    {
        size_t tot_read = 0, readSize = 0;
        writeHeader(outputFile, id, blockSize, isLast);
        while (tot_read < blockSize)
        {
            readSize = MIN(bufLen, blockSize - tot_read);
            handle_reading(buffer, readSize, inputFile);

            // Write to the output file
            safeWrite(buffer, 1, readSize, outputFile);
            tot_read += readSize;
        }
        safeFlush(outputFile);
        assert(tot_read == blockSize);

        // check signal
        bytesRead = read(fd, signal, strlen(TIMEOUTSIG));
        if (bytesRead > 0)
        {
            signal[bytesRead] = '\0'; // Null-terminate

            if (strcmp(signal, TIMEOUTSIG) == 0)
            {
                // log_with_timestamp("[ingate] Timeout SIG received. Processing timeout...");

                // keep reading until EOF
                readHeader(inputFile, &id, &blockSize, &isLast);
                while (!feof(inputFile))
                {
                    // printFileName(inputFile);
                    size_t tot_read = 0, readSize = 0;
                    writeHeader(outputFile, id, blockSize, isLast);
                    while (tot_read < blockSize)
                    {
                        readSize = MIN(bufLen, blockSize - tot_read);
                        handle_reading(buffer, readSize, inputFile);

                        // Write to the output file
                        safeWrite(buffer, 1, readSize, outputFile);
                        tot_read += readSize;
                    }
                    safeFlush(outputFile); // Ensure data is written to disk
                    assert(tot_read == blockSize);
                    readHeader(inputFile, &id, &blockSize, &isLast);
                }
                // close the inputFile since upstream EOF
                safeClose(inputFile);
                // ideally close the signal but cannot control if upstream has closed
                // close(fd)

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
            } // if signal == TIMEOUTSIG
        } // if bytesRead > 0

        readHeader(inputFile, &id, &blockSize, &isLast);
    }

    safeClose(inputFile);

    // clean up
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
    char* lambda_script_id = argv[4];

    inGate(inputFileBaseName, inputKeyBaseName, outputFileName, lambda_script_id);

    // clean up
    return 0;
}