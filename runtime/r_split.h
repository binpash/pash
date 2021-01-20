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

#define MINCHUNKS 12 //minimum number of chunks

#define  MIN(a,b) a > b ? b : a