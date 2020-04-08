#include <assert.h>
#include "eager_lib.h"

#define READ_WRITE_BUFFER_SIZE 2 * 1024

void EagerLoop(char* input, char* output, char* intermediate) {

    int doneReading = 0;
    int doneWriting = 0;

    struct timeval ts1, ts2, ts3, ts4;
    gettimeofday(&ts1, NULL);

    // Open the intermediate file for both reading and writing
    // TODO: Think of using O_TMPFILE
    int intermediateWriter = safeOpen3(intermediate, O_CREAT | O_WRONLY, S_IRWXU);
    int intermediateReader = safeOpen(intermediate, O_RDONLY);

    // It is fine for the input to block, since when we ask for it we
    // don't have anything in the intermediate file or buffer.
    int inputFd = safeOpen(input, O_RDONLY);
    debug("opened input file %s\n", input);

    debug("will open outputFile from %s \n", output);
    int outputFd = tryOpenOutput(output);
    while(outputFd < 0 && !doneReading) {
        if (readInputWriteToFile(inputFd, intermediateWriter, READ_WRITE_BUFFER_SIZE) == 0) {
            /* printf("Input was done before even output was opened\n"); */
            doneReading = 1;
        }
        outputFd = tryOpenOutput(output);
    }

    gettimeofday(&ts2, NULL);
    debug("Reading before the output was opened took %lu us\n",
           (ts2.tv_sec - ts1.tv_sec) * 1000000 + ts2.tv_usec - ts1.tv_usec);

    if (doneReading) {
        debug("will block to open outputFile from %s \n", output);
        outputFd = blockOpenOutput(output);
    }

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
            if (readInputWriteToFile(inputFd, intermediateWriter, READ_WRITE_BUFFER_SIZE) == 0) {
                doneReading = 1;
                debug("Input is done!\n");
                break;
            }
        }
        if (FD_ISSET(outputFd, &writeFds)) {
            // There is something to read in the intermediate file. We
            // then read from the intermediate file and write to the
            // output.
            off_t intermediateFileDiff =
                safeLseek(intermediateWriter) - safeLseek(intermediateReader);

            // There is nothing to read in the intermediate file. We
            // then have to block read from the input until it gives
            // us something.
            if (intermediateFileDiff == 0) {
                if (readInputWriteToFile(inputFd, intermediateWriter, READ_WRITE_BUFFER_SIZE) == 0) {
                    doneReading = 1;
                    debug("Input is done!\n");
                    break;
                }
                intermediateFileDiff =
                    safeLseek(intermediateWriter) - safeLseek(intermediateReader);
            }

            // Here we know that the intermediate file has something.
            assert(intermediateFileDiff > 0);
            safeWriteOutput(outputFd, intermediateReader, intermediateFileDiff, &doneWriting);
        }
    }

    gettimeofday(&ts3, NULL);
    debug("Select loop took %lu us\n",
           (ts3.tv_sec - ts2.tv_sec) * 1000000 + ts3.tv_usec - ts2.tv_usec);

    // If reading is done, close its file descriptor.
    if (doneReading) {
        close(inputFd);
    }

    // Set the output file to blocking so that writing to it doesn't lead to a busy loop
    fdSetBlocking(outputFd, 1);

    // Output the rest of the intermediate file
    outputRestIntermediateFile(outputFd, intermediateWriter, intermediateReader, &doneWriting);

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
