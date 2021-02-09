#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>

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
        PRINTDBG("Id read failed\n");
    }
    
    if((ret = fread(blockSize, sizeof(size_t), 1, inputFile)) < 0) {
        PRINTDBG("Blocksize read failed\n");
    };
}

void writeHeader(FILE* destFile, int64_t id, size_t blocksize) {
  fwrite(&id, sizeof(int64_t), 1, destFile);
  fwrite(&blocksize, 1, sizeof(size_t), destFile);
}