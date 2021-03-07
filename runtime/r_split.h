#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>

#ifdef DEBUG
#define PRINTDBG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define PRINTDBG(fmt, ...)
#endif

typedef __uint64_t uint64_t;

#define MAX_LINE_LENGTH 1e6

#define LOC __FILE__ // TODO expand

#define NUMOUTFILESLIMIT 4

#define CHUNKSIZE 1024*1024 //chunk will be bigger depending on where the line ends
#define BUFLEN BUFSIZ

#define MINCHUNKS 12 //minimum number of chunks

#define  MIN(a,b) a > b ? b : a

void readHeader(FILE* inputFile, int64_t *id, size_t *blockSize) {
    int ret;
    if ((ret = fread(id, sizeof(int64_t), 1, inputFile)) < 0) {
        fprintf(stderr, "Id read failed\n");
        exit(1);
    }
    
    if((ret = fread(blockSize, sizeof(size_t), 1, inputFile)) < 0) {
        fprintf(stderr, "Blocksize read failed\n");
        exit(1);
    };
}


void safeWrite(char* buffer, size_t bytes, size_t count, FILE* outputFile) {
    size_t len;
    if((len = fwrite(buffer, bytes, count, outputFile)) != count) {
      fprintf(stderr, "write failed count %ld, wrote %ld\n", count, len);
      exit(1);
    }
    fflush(outputFile);
}

void writeHeader(FILE* destFile, int64_t id, size_t blocksize) {
  safeWrite((char *) &id, sizeof(int64_t), 1, destFile);
  safeWrite((char *) &blocksize, sizeof(size_t), 1, destFile);
}