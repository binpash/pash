#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/select.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/sendfile.h>

/* /\* According to earlier standards *\/ */
/* #include <sys/time.h> */
/* #include <sys/types.h> */
/* #include <unistd.h> */

#ifndef __DEBUG__
#define __DEBUG__

#ifdef DEBUG
#define debug(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define debug(fmt, ...) ((void)0)
#endif

#endif

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define READ_WRITE_BUFFER_SIZE 1 * 1024

int safe_open3(const char *pathname, int flags, mode_t mode) {
    int fd = open(pathname, flags, mode);
    if (fd < 0) {
        printf("could not open file%s\n", pathname);
        exit(1);
    }
    return fd;
}

int safe_open(const char *pathname, int flags) {
    // The mode will be ignored anyway
    return safe_open3(pathname, flags, S_IRWXU);
}

int try_open_output(const char *pathname) {
    int outputFd = open(pathname, O_WRONLY | O_NONBLOCK);
    if (outputFd < 0) {
        if (errno == ENXIO) {
            // This means that noone has opened the output file for
            // reading, in that case we can read some of our input.
            /* printf("Noone has opened the output file for reading\n"); */
        } else {
            printf("could not open output file(s)\n");
            exit(1);
        }
    }
    return outputFd;
}

// Returns the number of bytes read, or 0 if the input was done.
int readInputWriteToFile(int inputFd, int intermediateWriter) {

    // TODO: Maybe allocate that buffer as a static or global to not
    // allocate it in the stack several times.
    ssize_t inputBytesRead = 0;
    ssize_t inputBytesWritten = 0;
    char inputBuf[READ_WRITE_BUFFER_SIZE];

    inputBytesRead = read(inputFd, inputBuf, sizeof(inputBuf));
    if (inputBytesRead < 0) {
        printf("Error: Couldn't read from input!\n");
        exit(1);
    }
    if (inputBytesRead == 0) {
        /* printf("Input is done!\n"); */
        return 0;
    }
    /* printf("Read %ld bytes from input\n", inputBytesRead); */

    inputBytesWritten = write(intermediateWriter, inputBuf, inputBytesRead);
    // TODO: I probably have to gracefully handle this case
    if (inputBytesWritten != inputBytesRead) {
        printf("Error: Didn't write all bytes to intermediate file!\n");
        exit(1);
    }
    return inputBytesRead;
}

// Returns the number of bytes written or 0 if the output is done
int writeOutput(int outputFd, const char* outputBuf, ssize_t bytesToWrite) {
    ssize_t bytesWritten =
        write(outputFd, outputBuf, bytesToWrite);
    if (bytesWritten < 0) {
        printf("Error: Couldn't write to output!\n");
        exit(1);
    }
    return bytesWritten;
}

// Returns 0 if output was done or a positive number otherwise
int emptyBuffer(int outputFd, const char* outputBuf, ssize_t* outputBytesRead, ssize_t* outputBytesWritten) {
    ssize_t newBytesWritten = 1;
    while(*outputBytesRead - *outputBytesWritten > 0) {
        newBytesWritten = writeOutput(outputFd, outputBuf, *outputBytesRead - *outputBytesWritten);
        if (newBytesWritten == 0) {
            debug("Output is done!\n");
            break;
        } else if (newBytesWritten < *outputBytesRead - *outputBytesWritten) {
            debug("didn't write everything\n");
        }
        *outputBytesWritten += newBytesWritten;
    }
    return newBytesWritten;
}

void outputRestIntermediateFile(int outputFd, int intermediateWriter, int intermediateReader,
                                char* outputBuf, int* doneWriting) {

    ssize_t outputBytesRead = 0;
    ssize_t outputBytesWritten = 0;
    // If writing is not done and there are things left in the buffer
    // or file, empty the buffer and intermediate files
    ssize_t intermediateFileBytesToOutput =
        lseek(intermediateWriter, 0, SEEK_CUR) - lseek(intermediateReader, 0, SEEK_CUR);
    // TODO: Is there a way to optimize this by just copying the rest
    // of the input in the output at once (e.g. by using cat)?
    while(!(*doneWriting) && intermediateFileBytesToOutput > 0) {

        // Fill the intermediate buffer
        if (intermediateFileBytesToOutput > 0) {
            outputBytesRead =
                read(intermediateReader, outputBuf,
                     MIN(intermediateFileBytesToOutput, sizeof(intermediateReader)));
            if (outputBytesRead < 0) {
                printf("Error: Didn't read from intermediate file!\n");
                exit(1);
            }
            outputBytesWritten = 0;
        }

        // Empty the intermediate buffer
        if (emptyBuffer(outputFd, outputBuf, &outputBytesRead, &outputBytesWritten) == 0) {
            *doneWriting = 1;
            break;
        }

        intermediateFileBytesToOutput =
            lseek(intermediateWriter, 0, SEEK_CUR) - lseek(intermediateReader, 0, SEEK_CUR);
    }

    return;
}

// WARNING: It is important to make sure that no operation (open, read,
// write).
void EagerLoop(char* input, char* output, char* intermediate) {

    int doneReading = 0;
    int doneWriting = 0;
    ssize_t outputBytesRead = 0;
    ssize_t outputBytesWritten = 0;
    char outputBuf[READ_WRITE_BUFFER_SIZE];

    struct timeval ts1, ts2, ts3, ts4;
    gettimeofday(&ts1, NULL);

    // Open the intermediate file for both reading and writing
    // TODO: Think of using O_TMPFILE
    int intermediateWriter = safe_open3(intermediate, O_CREAT | O_WRONLY, S_IRWXU);
    int intermediateReader = safe_open(intermediate, O_RDONLY);

    // It is fine for the input to block, since when we ask for it we
    // don't have anything in the intermediate file or buffer.
    int inputFd = safe_open(input, O_RDONLY);
    debug("opened input file %s\n", input);

    debug("will open outputFile from %s \n", output);
    int outputFd = try_open_output(output);
    while(outputFd < 0) {
        if (readInputWriteToFile(inputFd, intermediateWriter) == 0) {
            /* printf("Input was done before even output was opened\n"); */
            doneReading = 1;
        }
        outputFd = try_open_output(output);
    }

    gettimeofday(&ts2, NULL);
    debug("Reading before the output was opened took %lu us\n",
           (ts2.tv_sec - ts1.tv_sec) * 1000000 + ts2.tv_usec - ts1.tv_usec);

    // TODO: Optimize this loop by using sendfile
    while (!doneReading && !doneWriting) {
        fd_set readFds;
        fd_set writeFds;
        int maxFd;

        FD_ZERO(&readFds); // Clear FD set for select
        FD_ZERO(&writeFds); // Clear FD set for select
        FD_SET(inputFd, &readFds);
        FD_SET(outputFd, &writeFds);

        maxFd = MAX(inputFd, outputFd);

        // TODO: Should I handle some error here?
        select(maxFd + 1, &readFds, &writeFds, NULL, NULL);

        // TODO: Make writing more preferable by trying to rewrite if
        // it was possible once since this would help not accumulate
        // memory in this intermediate buffer.
        if (FD_ISSET(inputFd, &readFds)) {
            if (readInputWriteToFile(inputFd, intermediateWriter) == 0) {
                doneReading = 1;
                debug("Input is done!\n");
                break;
            }
        }
        if (FD_ISSET(outputFd, &writeFds)) {
            // If there is nothing in our intermediate buffer, we have to fill it.
            if (outputBytesRead == outputBytesWritten) {
                // 1. The intermediate buffer is empty but there is
                // something to read in the intermediate file. We then
                // read from the intermediate file and write to the
                // output.
                int intermediateFileDiff =
                    lseek(intermediateWriter, 0, SEEK_CUR) - lseek(intermediateReader, 0, SEEK_CUR);

                // 2. There is nothing to read in the intermediate
                // file. We then have to block read from the input until
                // it gives us something.
                if (intermediateFileDiff == 0) {
                    if (readInputWriteToFile(inputFd, intermediateWriter) == 0) {
                        doneReading = 1;
                        debug("Input is done!\n");
                        break;
                    }
                    intermediateFileDiff =
                        lseek(intermediateWriter, 0, SEEK_CUR) - lseek(intermediateReader, 0, SEEK_CUR);
                }

                // Here we know that the intermediate file has something.
                assert(intermediateFileDiff > 0);

                // Alternative 1
                outputBytesRead =
                    read(intermediateReader, outputBuf,
                         MIN(intermediateFileDiff, sizeof(outputBuf)));
                if (outputBytesRead <= 0) {
                    printf("Error: Didn't read from intermediate file!\n");
                    exit(1);
                }
                outputBytesWritten = 0;

                /* // Alternative 2 */
                /* ssize_t intermediateFileBytesToOutput; */
                /* ssize_t res; */
                /* intermediateFileBytesToOutput = */
                /*     lseek(intermediateWriter, 0, SEEK_CUR) - lseek(intermediateReader, 0, SEEK_CUR); */
                /* res = sendfile(outputFd, intermediateReader, 0, intermediateFileBytesToOutput); */
                /* if (res < 0 && errno != EAGAIN) { */
                /*     printf("ERROR: %s, when outputing!\n", strerror(errno)); */
                /*     exit(1); */
                /* } */
            }

            // Now the intermediate buffer certainly has data
            assert(outputBytesRead - outputBytesWritten > 0);
            ssize_t newBytesWritten =
                writeOutput(outputFd, outputBuf, outputBytesRead - outputBytesWritten);
            if (newBytesWritten == 0) {
                debug("We tried to write %ld, but output is done!\n", outputBytesRead - outputBytesWritten);
                doneWriting = 1;
                break;
            }
            outputBytesWritten += newBytesWritten;
        }
    }

    gettimeofday(&ts3, NULL);
    debug("Select loop took %lu us\n",
           (ts3.tv_sec - ts2.tv_sec) * 1000000 + ts3.tv_usec - ts2.tv_usec);

    // If reading is done, close its file descriptor.
    if (doneReading) {
        close(inputFd);
    }

    // Output the intermediate buffer
    if (emptyBuffer(outputFd, outputBuf, &outputBytesRead, &outputBytesWritten) == 0) {
        doneWriting = 1;
    }

    // Alternative 1: Output the rest of the intermediate file
    /* outputRestIntermediateFile(outputFd, intermediateWriter, intermediateReader, */
    /*                            outputBuf, &doneWriting); */

    // (Faster) Alternative 2: Output the rest of the file using sendfile
    // TODO: Put that in a function
    ssize_t intermediateFileBytesToOutput;
    ssize_t res;
    do {
        intermediateFileBytesToOutput =
            lseek(intermediateWriter, 0, SEEK_CUR) - lseek(intermediateReader, 0, SEEK_CUR);
        res = sendfile(outputFd, intermediateReader, 0, intermediateFileBytesToOutput);
        if (res < 0 && errno != EAGAIN) {
            printf("ERROR: %s, when outputing!\n", strerror(errno));
            exit(1);
        }
    } while (0 != res && res < intermediateFileBytesToOutput);

    if (res == 0) {
        // TODO: This might happen if output is a `head` or sth. I am
        // not sure what should we do in that case.
        doneWriting = 1;
    }

    gettimeofday(&ts4, NULL);
    debug("Finishing up writing the intermediate file took %lu us\n",
           (ts4.tv_sec - ts3.tv_sec) * 1000000 + ts4.tv_usec - ts3.tv_usec);


    // TODO: We have to handle the case where reading is not done but
    // writing is. In that case, we would like to draw all the input
    // to /dev/null (e.g. if our producer is tee).
    if(!doneReading) {
        printf("Error: Reading was not done when writing was!\n");
        exit(1);
    }
    close(outputFd);
    close(intermediateReader);
    close(intermediateWriter);
}

int main(int argc, char* argv[]) {
    // arg#1 -> input file name
    // arg#2 -> output file name
    // arg#3 -> intermediate file name
    if (argc < 4) {
        printf("ERROR: missing input!\n");
        exit(1);
    }

    char* inputFileName = calloc(strlen(argv[1]) + 1, sizeof(char));
    strcpy(inputFileName, argv[1]);

    char* outputFileName = calloc(strlen(argv[2]) + 1, sizeof(char));
    strcpy(outputFileName, argv[2]);
    char* intermediateFileName = calloc(strlen(argv[3]) + 1, sizeof(char));
    strcpy(intermediateFileName, argv[3]);

    EagerLoop(inputFileName, outputFileName, intermediateFileName);

    return 0;
}
