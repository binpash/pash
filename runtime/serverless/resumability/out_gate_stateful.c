#include "r_split.h"
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>

volatile int timer_expired = 0;
#define BUFFER_SIZE 100000

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
}

void outGate(FILE *inputFile, char* outputFileBaseName, 
            char* outputKeyBaseName, FILE *signalFile, 
            int timeout, char* lambda_script_id)
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
    size_t bufLen = BUFFER_SIZE;
    char *buffer = malloc(bufLen + 1);
    char *incompleteLine = malloc(bufLen + 1);
    ssize_t len = 0, headSize = 0, restSize = 0, prevRestSize = 0;
    int sentNewline = 1;

    while ((len = fread(buffer, 1, bufLen, inputFile)) > 0 )
    {
        // print it as [outgate outputkeybase]
        // printf("[outgate %s] Read %ld bytes\n", outputKeyBaseName, len);
        // find newline character
        for ( int i= len-1; i>=(len-1)/2; i--)
        {
            if (buffer[i] == '\n') 
            {
                headSize = i + 1;
                restSize = len - headSize;
                break;
            }
        }

        // no newline found
        if (headSize == 0) 
        {
            // send the incompleteLine and the current buffer
            headSize = len;
            if (prevRestSize)
                safeWrite(incompleteLine, 1, prevRestSize, outputFile);
            safeWriteWithFlush(buffer, 1, headSize, outputFile);
            // printf("[outgate %s] Send non-line: incomplete with size %ld, buffer with size %ld\n", outputKeyBaseName, prevRestSize, headSize);
            prevRestSize = 0;
            sentNewline = 0;
        } else {
            // send up to newline character, save the rest
            if (prevRestSize)
                safeWrite(incompleteLine, 1, prevRestSize, outputFile);
            safeWriteWithFlush(buffer, 1, headSize, outputFile);
            // printf("[outgate %s] Send line: incomplete with size %ld, buffer with size %ld\n", outputKeyBaseName, prevRestSize, headSize);
            memcpy(incompleteLine, buffer + headSize, restSize);
            prevRestSize = restSize;
            sentNewline = 1;
        }

        headSize = restSize = 0;

        // check the timer
        if (timer_expired)
        {
            // start redirecting if finish sending a line
            if (sentNewline)
            {
                // send signal to ingate
                safeWriteWithFlush(TIMEOUTSIG, sizeof(char), strlen(TIMEOUTSIG), signalFile);

                sleep(1); // wait for ingate to finish reading

                // close current output file
                safeClose(outputFile);

                // construct new tcp to new lambda
                version += 1;
                char* tcpKey = combineStrings(outputKeyBaseName, "_v", version);
                char* outputFileName = combineStrings(outputFileBaseName, "_v", version);
                createAndOpenFIFO(outputFileBaseName, "_v", version, &outputFile, 0666, O_RDWR, "w");

                // should i waitpid first? this might be slower cuz we just need to invoke the lambda faster to consume input
                // if (data_pid != -1)
                //     waitpid(data_pid, NULL, 0);
                pashSend(tcpKey, outputFileName);

                // invoke lambd
                invokeLambda(lambda_script_id, version);

                // clear timer, start a new timer for new lambda
                timer_expired = 0;
                start_timer(timeout);

            } else {
                // need to continue because newline is not sent
                continue;
            } 
        }
    }
    // printf("[outgate %s] Finished reading input file\n", outputKeyBaseName);

    if (len < 0)
    {
        perror(LOC);
        exit(1);
    }

    if (prevRestSize > 0)
    {
        // printf("[outgate %s] Sending incomplete line\n", outputKeyBaseName);
        safeWriteWithFlush(incompleteLine, 1, prevRestSize, outputFile);
    }

    
    // close output file
    safeClose(outputFile);

    // clean up
    free(buffer);
    free(incompleteLine);
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

    char* outputFileBaseName = argv[2];
    char* outputKeyBaseName = argv[3];

    // FILE *signalFile = fopen(argv[4], "w");
    // if (!signalFile)
    // {
    //     perror("Error opening signal file");
    //     exit(1);
    // }
    
    int timeout = atoi(argv[4]);

    char* lambda_script_id = argv[5];

    FILE *signalFile = NULL;

    createAndOpenFIFO(lambda_script_id, "_sig_", 0, &signalFile, 0666, O_WRONLY, "w");

    // log_with_timestamp("[outgate %s] Start processing data with timeout %d", outputKeyBaseName, timeout);

    outGate(inputFile, outputFileBaseName, outputKeyBaseName,signalFile, timeout, lambda_script_id);

    // close files
    safeClose(inputFile);
    safeClose(signalFile);
    // remove signal file
    char* signalFileName = combineStrings(lambda_script_id, "_sig_", 0);
    remove(signalFileName);

    return 0;
}