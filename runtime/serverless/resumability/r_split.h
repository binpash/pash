#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <err.h>

#ifdef DEBUG
#define PRINTDBG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define PRINTDBG(fmt, ...)
#endif

typedef __uint64_t uint64_t;
typedef int8_t bool;
#define MAX_LINE_LENGTH 1e6

#define LOC __FILE__ // TODO expand

#define NUMOUTFILESLIMIT 4

#define CHUNKSIZE 1024 * 1024 //chunk will be bigger depending on where the line ends
#define BUFLEN 4*BUFSIZ

#define MINCHUNKS 12 //minimum number of chunks

#define MIN(a, b) a > b ? b : a
#define MAX(a, b) a < b ? b : a

typedef struct block_header {
  int64_t id;
  size_t blockSize;
  int8_t isLast; 
} block_header;


void safeWriteWithFlush(char *buffer, size_t bytes, size_t count, FILE *outputFile)
{
  size_t len;
  if ((len = fwrite(buffer, bytes, count, outputFile)) != count)
  {
    // Reader terminated early 
    // Todo: might need to ignore sigpipe signal
    if (len == -1 && errno == EPIPE) {
      exit(EXIT_SUCCESS);
    }
    err(2, "write failed count %lu, wrote %lu", count, len);
  }
  fflush(outputFile);
}
void safeWrite(char *buffer, size_t bytes, size_t count, FILE *outputFile)
{
  size_t len;
  if ((len = fwrite(buffer, bytes, count, outputFile)) != count)
  {
    // Reader terminated early
    if (len == -1 && errno == EPIPE) {
      exit(EXIT_SUCCESS);
    }

    err(2, "write failed count %lu, wrote %lu", count, len);
  }
}

// Clear block in case of broken pipe
void clear_block(FILE *src, size_t tot_read, size_t blockSize)
{
  char* buffer = malloc(BUFLEN);
  while (tot_read < blockSize)
  {
    size_t readSize = MIN(BUFLEN, blockSize - tot_read);
    // fprintf(stderr, "reading stdin\n");
    if (readSize && fread(buffer, 1, readSize, src) != readSize)
    {
      err(2, "There is a problem with reading the block");
    }
  }
  free(buffer);
}

//taken from https://github.com/dspinellis/dgsh/blob/master/core-tools/src/dgsh-tee.c
void non_block_fd(int fd, const char *name)
{
  int flags = fcntl(fd, F_GETFL, 0);
  if (flags < 0)
    err(2, "Error getting flags for %s", name);
  if (fcntl(fd, F_SETFL, flags | O_NONBLOCK) < 0)
    err(2, "Error setting %s to non-blocking mode", name);
}

void block_fd(int fd, const char *name)
{
  int flags = fcntl(fd, F_GETFL, 0);
  if (flags < 0)
    err(2, "Error getting flags for %s", name);
  if (fcntl(fd, F_SETFL, flags & (~O_NONBLOCK)) < 0)
    err(2, "Error setting %s to blocking mode", name);
}

void readHeader(FILE *inputFile, int64_t *id, size_t *blockSize, bool *isLast)
{
  size_t ret;
  block_header header;
  if ((ret = fread(&header, sizeof(block_header), 1, inputFile)) < 0)
    err(2, "Reading header failed");
  *id = header.id;
  *blockSize = header.blockSize;
  *isLast = header.isLast;
}

void writeHeader(FILE *destFile, int64_t id, size_t blockSize, bool isLast)
{
  block_header header = {id, blockSize, isLast};
  safeWrite((char *)&header, sizeof(block_header), 1, destFile);
}

size_t handle_reading(char *buffer, size_t size, FILE* stream) {
  size_t len = 0;
  // size_t readSize = MIN(bufLen, blockSize-total_read);
  if (size && (len = fread(buffer, 1, size, stream)) != size)
  {   
    // This fixes a problem that happens sometimes when pash kill signal effects the process we were reading from but not the current process
    // Might need to extend this to gracefully exit if errno == SIGPIPE
    if (feof(stream) && !errno) {
        PRINTDBG("Pipe closed before the full block data was written");
        exit(0);
    } else {
        err(2, "There is a problem with reading the block");
    }
  }
  return len;
}
