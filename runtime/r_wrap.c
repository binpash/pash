#include "r_split.h"
#include <unistd.h>
#include<sys/wait.h> 

#define READ_END 0
#define WRITE_END 1

void processCmd(char* args[]) {
    size_t bufLen = BUFLEN; //buffer length, would be resized as needed
    int64_t id;
    size_t blockSize;
    char* buffer = malloc(bufLen + 1);

    readHeader(stdin, &id, &blockSize);
    while(!feof(stdin)) {     
        int fdIn[2];
        int fdOut[2];
        if (pipe(fdIn) < 0) {
            perror("pipe failed");
        }

        if (pipe(fdOut) < 0) {
            perror("pipe failed");
        }
        
        int pid = fork();
        if (pid == 0) {
            dup2(fdOut[WRITE_END], STDOUT_FILENO);
            dup2(fdIn[READ_END], STDIN_FILENO);
            close(fdOut[READ_END]);
            close(fdOut[WRITE_END]);
            close(fdIn[READ_END]);
            close(fdIn[WRITE_END]);
            free(buffer);
            execvp(args[0], args);
            //shouldn't get here
            perror("Exec failed");
        } else {
            close(fdIn[READ_END]);
            close(fdOut[WRITE_END]);
            
            FILE* execOutFile = fdopen(fdOut[READ_END], "r");
            //Read batch
            size_t tot_read = 0, readSize = 0;
            while (tot_read < blockSize) {
                readSize = MIN(bufLen, blockSize-tot_read);
                if (fread(buffer, 1, readSize, stdin) != readSize) {
                    fprintf(stderr, "r_wrap: There is a problem with reading the block\n");
                    exit(1);
                }
                
                //Write to forked process
                if (write(fdIn[WRITE_END], buffer, readSize) != readSize) {
                    fprintf(stderr, "Failed wrtiting to pipe\n");
                    exit(1);
                }

                tot_read += readSize;
            }
            close(fdIn[WRITE_END]);
            
            assert(tot_read == blockSize);

            //read output of forked process (do I need to wait or is read blocking enough?)
            waitpid(pid, NULL, 0);
            size_t outBufLen = 0, len = 0, currLen = 0;
            char *cmdOutput = NULL;
            while ((len = fread(buffer, 1, bufLen, execOutFile)) > 0) {
                if ((currLen + len) > outBufLen) {
                    outBufLen = currLen + len + CHUNKSIZE;
                    cmdOutput = realloc(cmdOutput, outBufLen + 1);
                }
                memcpy(cmdOutput + currLen, buffer, len);
                currLen += len;
            }
            fclose(execOutFile);

            //write block to stdout
            writeHeader(stdout, id, currLen);
            if (fwrite(cmdOutput, 1, currLen, stdout) != currLen) {
                perror("fwrite failed");
            }

            //update header (ordered at the end so !feof works) and cleanup
            readHeader(stdin, &id, &blockSize);
            free(cmdOutput);
        }  
        
    }
    free(buffer);
}

int main(int argc, char* argv[]) {
    //arg1: command
    //args 2.. : arguments for command
    //input is from stdin, out to stdout
    char **args = NULL;
    if (argc < 2) {
		/* default behavior is to echo all the filenames */
		fprintf(stderr, "missing input!\n");
        exit(1);
    }

    //process arguments
    args = malloc(sizeof(char *)*(argc-1));
    for (int i = 1; i < argc; i++) {
       args[i-1] = malloc(strlen(argv[i]));
       strcpy(args[i-1], argv[i]);
    }

    processCmd(args);
}
