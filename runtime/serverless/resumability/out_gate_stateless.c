#include "r_split.h"
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>

volatile int timer_expired = 0;

void timer_handler(int signo)
{
    if (signo == SIGALRM)
    {
        timer_expired = 1;
    }
}

void start_timer(const int timeout)
{
    struct itimerval timer;
    timer.it_value.tv_sec = timeout;
    timer.it_value.tv_usec = 0;
    timer.it_interval.tv_sec = timeout;
    timer.it_interval.tv_usec = 0;

    if (setitimer(ITIMER_REAL, &timer, NULL) < 0)
    {
        perror("Failed to set timer");
        exit(1);
    }
    // printf("[outgate] Timer started with interval: %d seconds.\n", timeout);
}

void outGate(FILE *inputFile, char *outputFileBaseName,
             char *outputKeyBaseName, FILE *signalFile,
             int timeout, char *lambda_script_id)
{
    // start timer
    signal(SIGALRM, timer_handler);
    start_timer(timeout);

    // construct output file name
    int version = 0;
    char *outputFileName = combineStrings(outputFileBaseName, "_v", version);
    FILE *outputFile = fopen(outputFileName, "w");
    if (!outputFile)
    {
        perror("Error opening output file");
        exit(1);
    }

    // start processing data
    size_t bufLen = BUFLEN;
    char *buffer = malloc(bufLen + 1);
    int64_t id;
    size_t blockSize;
    bool isLast;

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

        // check the timer
        if (timer_expired)
        {
            // send signal to ingate
            safeWriteWithFlush(TIMEOUTSIG, sizeof(char), strlen(TIMEOUTSIG), signalFile);

            // close current output file
            safeClose(outputFile);

            // construct new tcp to new lambda
            version += 1;
            char *tcpKey = combineStrings(outputKeyBaseName, "_v", version);
            char *outputFileName = combineStrings(outputFileBaseName, "_v", version);
            createAndOpenFIFO(outputFileBaseName, "_v", version, &outputFile, 0666, O_RDWR, "w");

            // should i waitpid first? this might be slower cuz we just need to invoke the lambda faster to consume input
            // if (data_pid != -1)
            //     waitpid(data_pid, NULL, 0);
            pashSend(tcpKey, outputFileName);


            invokeLambda(lambda_script_id, version);
            // clear timer, start a new timer for new lambda
            timer_expired = 0;
            start_timer(timeout);
        }
        readHeader(inputFile, &id, &blockSize, &isLast);
        // printf("[outgate] Reading block %ld size %ld | eof %d \n", id, blockSize, feof(inputFile));
    }

    // clean up
    free(buffer);
}

// arg#1: data input file
// arg#2: data output file base name
// arg#3: data output key base name
// arg#4: output signal file
// arg#5: timeout
int main(int argc, char *argv[])
{
    FILE *inputFile = fopen(argv[1], "r");
    if (!inputFile)
    {
        perror("Error opening input file");
        exit(1);
    }

    char *outputFileBaseName = argv[2];
    char *outputKeyBaseName = argv[3];

    int timeout = atoi(argv[4]);
    char* lambda_script_id = argv[5];

    FILE *signalFile = NULL;
    
    createAndOpenFIFO(lambda_script_id, "_sig_", 0, &signalFile, 0666, O_WRONLY, "w");

    outGate(inputFile, outputFileBaseName, outputKeyBaseName, signalFile, timeout, lambda_script_id);

    // close files
    safeClose(inputFile);
    safeClose(signalFile);
    // remove the signal file
    char *signalFileName = combineStrings(lambda_script_id, "_sig_", 0);
    remove(signalFileName);

    return 0;
}