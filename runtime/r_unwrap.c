#include "r_split.h"

void unwrap(FILE* inputFile) {
    size_t bufLen = BUFLEN; //buffer length, would be resized as needed
    int64_t id;
    size_t blockSize;
    char* buffer = malloc(bufLen + 1);

    readHeader(stdin, &id, &blockSize);
    while(!feof(stdin)) {     
        //Read batch
        size_t tot_read = 0, readSize = 0;
        while (tot_read < blockSize) {
            readSize = MIN(bufLen, blockSize-tot_read);
            if (fread(buffer, 1, readSize, stdin) != readSize) {
                fprintf(stderr, "r_wrap: There is a problem with reading the block\n");
                exit(1);
            }
            
            //Write to forked process
            safeWriteWithFlush(buffer, 1, readSize, stdout);
            tot_read += readSize;
        }
        // fflush(stdout);
        assert(tot_read == blockSize);

        //update header (ordered at the end so !feof works) and cleanup
        readHeader(stdin, &id, &blockSize);
    }  
        
    free(buffer);   
}

int main(int argc, char* argv[]) {
    //no arguments needed, can accept one argument for input file instead of stdin
    //input is from stdin, out to stdout
    FILE *inputFile = stdin; //defualt is stdin
    if (argc > 1) {
        inputFile = fopen(argv[1], "r");
    }

    unwrap(inputFile);

    //cleanup
    fclose(inputFile);
}
