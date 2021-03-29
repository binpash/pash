#include "r_split.h"
#include <unistd.h>
#include <sys/wait.h>

#define READ_END 0
#define WRITE_END 1

/*
Batch sizes in and out should be around the same size or smaller ideally. 
The code was not tested in cases were the fork would output significantly more lines 
than the input but that would definitly lead to incread memory cosumption in the whole pipeline.
*/
void processCmd(char *args[])
{
    size_t bufLen = BUFLEN, stdoutBlockBufLen = BUFLEN, writeBufLen = 2*BUFLEN; //buffer length, would be resized as needed
    int64_t id;
    size_t blockSize;
    char *buffer = malloc(bufLen + 1);
    char *readBuffer = malloc(bufLen + 1);
    char *writebuffer = malloc(writeBufLen + 1);
    char *stdoutBlock = malloc(stdoutBlockBufLen+1);
    //select
    fd_set readFds;
    fd_set writeFds;
    int maxFd;

    readHeader(stdin, &id, &blockSize);
    while (!feof(stdin))
    {
        int fdIn[2];
        int fdOut[2];
        if (pipe(fdIn) < 0)
        {
            perror("pipe failed");
            exit(1);
        }

        if (pipe(fdOut) < 0)
        {
            perror("pipe failed");
            exit(1);
        }

        int pid = fork();
        if (pid == 0)
        {
            dup2(fdOut[WRITE_END], STDOUT_FILENO);
            dup2(fdIn[READ_END], STDIN_FILENO);

            close(fdOut[READ_END]);
            close(fdOut[WRITE_END]);
            close(fdIn[READ_END]);
            close(fdIn[WRITE_END]);
            // free(buffer); copy-on-write fork optimize this better(freeing is more harm than good here)
            execvp(args[0], args);
            //shouldn't get here
            perror("Exec failed");
            exit(1);
        }
        else
        {
            close(fdIn[READ_END]);
            close(fdOut[WRITE_END]);

            int inputFd = fdOut[READ_END];
            int outputFd = fdIn[WRITE_END];

            // FILE *execOutFile = fdopen(inputFd, "rb");
            // FILE *execInFile = fdopen(outputFd, "wb");
            
            non_block_fd(inputFd, "r-wrap-fork-read");
            non_block_fd(outputFd, "r-wrap-fork-write");

            size_t len = 0, currReadLen = 0, currWriteLen = 0;

            //select prep
            maxFd = MAX(inputFd, outputFd);

            //Read batch
            size_t tot_read = 0, readSize = 0;
            while (tot_read < blockSize || currWriteLen > 0)
            {
                readSize = MIN(bufLen, blockSize - tot_read);
                // fprintf(stderr, "reading stdin\n");
                if (readSize && fread(buffer, 1, readSize, stdin) != readSize)
                {
                    err(2, "r_wrap: There is a problem with reading the block");
                }
                if ((currWriteLen + readSize) > writeBufLen) {
                    writeBufLen = 2*(currWriteLen + readSize);
                    writebuffer = realloc(writebuffer, writeBufLen+1);
                }
                memcpy(writebuffer + currWriteLen, buffer, readSize);
                currWriteLen += readSize;

                FD_ZERO(&writeFds); // Clear FD set for select
                while (!FD_ISSET(outputFd, &writeFds))
                {
                    // fprintf(stderr, "looping select\n");
                    FD_ZERO(&readFds);  // Clear FD set for select
                    FD_ZERO(&writeFds); // Clear FD set for select
                    FD_SET(inputFd, &readFds);
                    FD_SET(outputFd, &writeFds);

                    // Theoretically we don't need select because both pipes are blocking 
                    // but it saves cpu cycles
                    if (select(maxFd + 1, &readFds, &writeFds, NULL, NULL) < 0) {
                        if (errno == EINTR)
				            continue;
			            perror("select");
                    }
                    
                    if (FD_ISSET(inputFd, &readFds))
                    {
                        //Try reading from forked processs, nonblocking
                        if ((len = read(inputFd, readBuffer, bufLen)) > 0) {
                            if ((currReadLen + len) > stdoutBlockBufLen)
                            {
                                stdoutBlockBufLen = 2*(currReadLen + len);
                                stdoutBlock = realloc(stdoutBlock, stdoutBlockBufLen + 1);
                            }
                            memcpy(stdoutBlock + currReadLen, readBuffer, len);
                            currReadLen += len;
                        } else if (len < 0) {
                            switch (errno) {
					        case EAGAIN:
                                break;
                            default:
                                err(2, "r_wrap: failed reading from fork, error %d", errno);
                            }
                        }
                    }
                }
                // fprintf(stderr, "writing %ld bytes\n", readSize);
                // Write to forked process
                // For some reason nonblocking fd is not working correctly with fwrite. 
                // Check shortest-scripts.sh with big batch size to recreate the problem
                if ((len = write(outputFd, writebuffer, currWriteLen)) > 0) {
                    currWriteLen -= len; 
                    memmove(writebuffer, writebuffer + len, currWriteLen);
                } else if (len < 0) {
                    switch (errno) {
                        
					case EAGAIN:
						len = 0;
						break;
                    // TODO handle broken pipe if needed
                    case EPIPE:
					default:
						err(2, "r_wrap: error writing to fork, error %d", errno);
					}            
                }

                tot_read += readSize;
            }
            close(outputFd);
            // fprintf(stderr, "finished one fork\n");
            assert(tot_read == blockSize);

            // read output of forked process
            // block to make sure you read until the process exits and to save cpu cycles
            block_fd(inputFd, "r-wrap-fork-read");
            bool reached_eof = 0;
            while (!reached_eof) {
                if ((len = read(inputFd, buffer, bufLen)) > 0) {
                    if ((currReadLen + len) > stdoutBlockBufLen)
                    {
                        stdoutBlockBufLen = currReadLen + len + CHUNKSIZE;
                        stdoutBlock = realloc(stdoutBlock, stdoutBlockBufLen + 1);
                    }
                    memcpy(stdoutBlock + currReadLen, buffer, len);
                    currReadLen += len;
                } else {
                    if (len == 0)
                        reached_eof = 1;
                    else {
                        switch (errno) {
                            case EAGAIN:
                                break;
                            default:
                                err(2, "r_wrap: error writing to fork, error %d", errno);
                        }
                    }
                }
            }
            close(inputFd);

            //write block to stdout
            writeHeader(stdout, id, currReadLen);
            safeWriteWithFlush(stdoutBlock, 1, currReadLen, stdout);

            //update header (ordered at the end so !feof works) and cleanup
            readHeader(stdin, &id, &blockSize);
        }
        kill(pid, SIGKILL);
    }
    free(buffer);
    free(readBuffer);
    free(stdoutBlock);
    free(writebuffer);
}

int main(int argc, char *argv[])
{
    //arg1: command
    //args 2.. : arguments for command
    //input is from stdin, out to stdout
    char **args = NULL;
    if (argc < 2)
    {
        /* default behavior is to echo all the filenames */
        fprintf(stderr, "missing input!\n");
        exit(1);
    }

    //process arguments
    args = malloc(sizeof(char *) * (argc));
    for (int i = 1; i < argc; i++)
    {
        args[i - 1] = malloc(strlen(argv[i]) + 1);
        strcpy(args[i - 1], argv[i]);
    }
    args[argc - 1] = '\0';
    processCmd(args);
}
