#define _GNU_SOURCE
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
#ifdef __linux__
#include <sys/sendfile.h>
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/uio.h>
#endif

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


int safeOpen3(const char *pathname, int flags, mode_t mode);

int safeOpen(const char *pathname, int flags);

off_t safeLseek(int fd);

void fdSetBlocking(int fd, int blocking);

int tryOpenOutput(const char *pathname);

int blockOpenOutput(const char *pathname);

int readInputWriteToFile(int inputFd, int intermediateWriter, int bufferSize);

/* int bufferedReadInputWriteToFile(int inputFd, int intermediateWriter, int bufferSize); */

/* int writeOutput(int outputFd, const char* outputBuf, ssize_t bytesToWrite); */

/* int emptyBuffer(int outputFd, const char* outputBuf, ssize_t* outputBytesRead, ssize_t* outputBytesWritten); */

/* void bufferedOutputRestIntermediateFile(int outputFd, int intermediateWriter, int intermediateReader, */
/*                                         char* outputBuf, int* doneWriting); */

ssize_t safeWriteOutput(int outputFd, int intermediateReader,
                        off_t intermediateFileDiff, int* doneWriting);

void outputRestIntermediateFile(int outputFd, int intermediateWriter,
                                int intermediateReader, int* doneWriting);
