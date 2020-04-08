#include "eager_lib.h"

int safeOpen3(const char *pathname, int flags, mode_t mode) {
    int fd = open(pathname, flags, mode);
    if (fd < 0) {
        printf("could not open file%s\n", pathname);
        exit(1);
    }
    return fd;
}

int safeOpen(const char *pathname, int flags) {
    // The mode will be ignored anyway
    return safeOpen3(pathname, flags, S_IRWXU);
}

off_t safeLseek(int fd) {
    off_t offset = lseek(fd, 0, SEEK_CUR);
    if (offset < 0) {
        printf("ERROR: %s, when lseek-ing!\n", strerror(errno));
        exit(1);
    }
    return offset;
}

/**
 * Copied from: http://code.activestate.com/recipes/577384-setting-a-file-descriptor-to-blocking-or-non-block/
 * Set a file descriptor to blocking or non-blocking mode.
 *
 * @param fd The file descriptor
 * @param blocking 0:non-blocking mode, 1:blocking mode
 *
 **/
void fdSetBlocking(int fd, int blocking) {
    /* Save the current flags */
    int flags = fcntl(fd, F_GETFL, 0);
    if (flags == -1) {
        printf("could not get file descriptor %d flags\n", fd);
        exit(1);
    }

    if (blocking) {
        flags &= ~O_NONBLOCK;
    } else {
        flags |= O_NONBLOCK;
    }
    int res = fcntl(fd, F_SETFL, flags);
    if (res == -1) {
        printf("could not set file descriptor %d to blocking\n", fd);
        exit(1);
    }
    return;
}

int tryOpenOutput(const char *pathname) {
    int outputFd = open(pathname, O_WRONLY | O_NONBLOCK);
    if (outputFd < 0) {
        // ENXIO means that noone has opened the output file for
        // reading, in that case we can read some of our input.
        if (errno != ENXIO) {
            printf("could not open output file(s)\n");
            exit(1);
        }
    }
    return outputFd;
}

int blockOpenOutput(const char *pathname) {
    int outputFd = open(pathname, O_WRONLY);
    if (outputFd < 0) {
        printf("could not open output file: %s\n", pathname);
        exit(1);
    }
    fdSetBlocking(outputFd, 0);
    return outputFd;
}

// Returns the number of bytes read, or 0 if the input was done.
int readInputWriteToFile(int inputFd, int intermediateWriter, int bufferSize) {

    ssize_t res = splice(inputFd, 0, intermediateWriter, 0, bufferSize, 0);
    if (res < 0) {
        printf("Error: Couldn't read from input!\n");
        exit(1);
    }
    return res;
}


int bufferedReadInputWriteToFile(int inputFd, int intermediateWriter, int bufferSize) {
    // TODO: Maybe allocate that buffer as a static or global to not
    // allocate it in the stack several times.
    ssize_t inputBytesRead = 0;
    ssize_t inputBytesWritten = 0;

    printf("Is it valid to allocate an array of variable size?\n");
    exit(1);
    char inputBuf[bufferSize];

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

void bufferedOutputRestIntermediateFile(int outputFd, int intermediateWriter, int intermediateReader,
                                        char* outputBuf, int* doneWriting) {

    ssize_t outputBytesRead = 0;
    ssize_t outputBytesWritten = 0;
    // If writing is not done and there are things left in the buffer
    // or file, empty the buffer and intermediate files
    ssize_t intermediateFileBytesToOutput =
        safeLseek(intermediateWriter) - safeLseek(intermediateReader);
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
            safeLseek(intermediateWriter) - safeLseek(intermediateReader);
    }

    return;
}

ssize_t safeWriteOutput(int outputFd, int intermediateReader,
                        int intermediateFileDiff, int* doneWriting) {
    ssize_t res;
    res = sendfile(outputFd, intermediateReader, 0, intermediateFileDiff);
    if (res < 0 && errno != EAGAIN) {
        printf("ERROR: %s, when outputing!\n", strerror(errno));
        exit(1);
    } else if (res == 0) {
        debug("We tried to write %d, but output is done!\n", intermediateFileDiff);
        *doneWriting = 1;
    }
    return res;
}

void outputRestIntermediateFile(int outputFd, int intermediateWriter,
                                int intermediateReader, int* doneWriting) {
    ssize_t finalOffset = safeLseek(intermediateWriter);
    ssize_t intermediateFileBytesToOutput;
    ssize_t res;
    do {
        intermediateFileBytesToOutput =
             finalOffset - safeLseek(intermediateReader);
        res = safeWriteOutput(outputFd, intermediateReader, intermediateFileBytesToOutput, doneWriting);
    } while (!(*doneWriting) && res < intermediateFileBytesToOutput);

    return;
}
