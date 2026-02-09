#include "r_split.h"
#include "stdbool.h"

void SplitByBytes(FILE *inputFile, int batchSize, FILE *outputFiles[], unsigned int numOutputFiles)
{
  int current_file_id = 0;
  int64_t id = 0;
  size_t len = 0;
  FILE *outputFile = outputFiles[current_file_id];

  char *buffer = malloc(batchSize + 1);

  // Do round robin copying of the input file to the output files
  // Each block has a header of "ID blockSize\n"
  while ((len = fread(buffer, 1, batchSize, inputFile)) > 0)
  {
    //write header
    writeHeader(outputFile, id, len, 1);

    //write blocks
    fwrite(buffer, 1, len, outputFile);

    current_file_id = (current_file_id + 1) % numOutputFiles;
    outputFile = outputFiles[current_file_id];
    id += 1;
  }

  if (len < 0)
  {
    perror(LOC);
    exit(1);
  }

  //clean up
  free(buffer);
}

int find_new_line_pivot(char *buffer, int start_pos, int end_pos, bool backward) {
  if (backward) {
    for (int i = end_pos; i >= start_pos; i--) {
      if (buffer[i] == '\n')
        return i;
    }
  }
  else {
    for (int i = start_pos; i <= end_pos; i++) {
      if (buffer[i] == '\n')
        return i;
    }
  }
  return -1;
}

void SplitByLines(FILE *inputFile, int batchSize, FILE *outputFiles[], unsigned int numOutputFiles, bool add_header)
{
  int current_file_id = 0;
  int64_t id = 0;
  size_t len = 0, headSize = 0, restSize = 0, prevRestSize = 0, blockSize = 0;
  FILE *outputFile = outputFiles[current_file_id];

  char *buffer = malloc(batchSize + 1);
  char *incompleteLine = malloc(batchSize + 1);


  // First read an initial batch of W * batchSize to make sure we split equally incase data size is small
  size_t full_payload = batchSize * numOutputFiles;
  char *init_buffer = malloc(full_payload + 1);
  if((len = fread(init_buffer, 1, full_payload, inputFile)) > 0) {
    int start_pos = 0;
    bool is_last = true;
    size_t init_batch_size = batchSize;
    // Total file size is small enough that we need a smaller batchsize
    if (len < full_payload) {
      init_batch_size = len/numOutputFiles;
    }

    for (current_file_id = 0; current_file_id < numOutputFiles; current_file_id++) {
      outputFile = outputFiles[current_file_id];
      int next_start = 0;
      if (current_file_id < numOutputFiles - 1) {
        // Process output for the first n - 1 nodes
        int end_pos = (current_file_id + 1) * init_batch_size;
        int pivot = find_new_line_pivot(init_buffer, start_pos, end_pos, true);
        // If no newline in the range [start_pos, end_pos] then blocksize would be 0
        if (pivot == -1) {
          next_start = start_pos;
          blockSize = 0;
        } else {
          next_start = pivot + 1;
          blockSize = next_start - start_pos;
        }

      } else {
        // Process output for last node
        blockSize = len - start_pos;
        if (init_buffer[len - 1] == '\n' || feof(inputFile)) {
          is_last = true;
        } else {
          is_last = false;
        }

      }

      if (add_header)
        writeHeader(outputFile, id, blockSize, is_last);

      safeWriteWithFlush(init_buffer + start_pos, 1, blockSize, outputFile);

      start_pos = next_start;

      if (is_last) {
        id += 1;
      }
    }

    // This will wrap around if the last chunk was complete
    // otherwise keep pointer on the last file
    if (is_last) {
      current_file_id = 0;
    } else {
      current_file_id = numOutputFiles - 1;
    }
  }
  free(init_buffer);

  // Do round robin copying of the input file to the output files
  // Each block has a header of "ID blockSize\n"
  while ((len = fread(buffer, 1, batchSize, inputFile)) > 0)
  {
    outputFile = outputFiles[current_file_id];

    //find pivot point for head and rest
    for (int i = len - 1; i >= (len - 1)/2; i--) //only search to the middle
    {
      if (buffer[i] == '\n')
      {
        headSize = i + 1;
        restSize = len - headSize;
        break;
      }
    }

    //no new line character
    if (headSize == 0)
    {
      headSize = len;
      blockSize = prevRestSize + headSize;
      if (add_header) {
          if (feof(inputFile))
            writeHeader(outputFile, id, blockSize, 1);
          else
            writeHeader(outputFile, id, blockSize, 0);
      }
      if (prevRestSize)
        safeWrite(incompleteLine, 1, prevRestSize, outputFile);
      safeWriteWithFlush(buffer, 1, headSize, outputFile);

      // Prepare next iteration
      prevRestSize = 0;
    }
    else
    {
      blockSize = prevRestSize + headSize;
      //write header
      if (add_header)
        writeHeader(outputFile, id, blockSize, 1);
      //write blocks
      if (prevRestSize)
        safeWrite(incompleteLine, 1, prevRestSize, outputFile);
      safeWriteWithFlush(buffer, 1, headSize, outputFile);
      //update incompleteLine to the current block
      memcpy(incompleteLine, buffer + headSize, restSize);

      // Prepare next iteration
      current_file_id = (current_file_id + 1) % numOutputFiles;
      prevRestSize = restSize;
      id += 1;
    }

    headSize = restSize = 0;
  }

  if (len < 0)
  {
    perror(LOC);
    exit(1);
  }

  if (prevRestSize > 0)
  {
    if (add_header)
      writeHeader(outputFile, id, prevRestSize, 1);
    safeWriteWithFlush(incompleteLine, 1, prevRestSize, outputFile);
  }

  //clean up
  free(buffer);
  free(incompleteLine);
}

void SplitByLinesRaw(FILE *inputFile, int batchSize, FILE *outputFiles[], unsigned int numOutputFiles)
{
  int current_file_id = 0, len = 0;
  size_t bufLen = 0;
  FILE* outputFile = outputFiles[current_file_id];

  char* buffer = NULL;

  // Do round robin copying of the input file to the output files without any headers
  while ((len = getline(&buffer, &bufLen, inputFile)) > 0) {
      safeWrite(buffer, 1, len, outputFile);
      current_file_id = (current_file_id + 1) % numOutputFiles;
      outputFile = outputFiles[current_file_id];
  }
  //clean up
  free(buffer);

}



void SplitInput(char *input, int batchSize, char *outputFileNames[], unsigned int numOutputFiles, bool useBytes, bool raw)
{
  PRINTDBG("%s: will split input\n", __func__);
  //TODO: find better way?
  FILE **outputFiles = malloc(sizeof(FILE *) * numOutputFiles);
  for (int i = 0; i < numOutputFiles; i++)
  {
    outputFiles[i] = fopen(outputFileNames[i], "w");
    if (!outputFiles[i])
    {
      perror(LOC);
      exit(1);
    }
  }

  FILE *inputFile = fopen(input, "r");
  if (!inputFile)
  {
    perror(LOC);
    exit(1);
  }
  PRINTDBG("%s: Opened input file %s\n", __func__, input);

  if (raw)
  {
    // SplitByLinesRaw(inputFile, batchSize, outputFiles, numOutputFiles);
    SplitByLines(inputFile, batchSize, outputFiles, numOutputFiles, 0);
  }
  else if (useBytes)
  {
    //if batchSize isn't set we can approximate it
    if (batchSize == 0)
    {
      int inputfd = fileno(inputFile);
      struct stat buf;
      fstat(inputfd, &buf);
      size_t inputSize = buf.st_size;
      batchSize = MIN(inputSize / MINCHUNKS, CHUNKSIZE); //autotune this better?
    }
    SplitByBytes(inputFile, batchSize, outputFiles, numOutputFiles);
  }
  else
  {
    SplitByLines(inputFile, batchSize, outputFiles, numOutputFiles, 1);
  }

  PRINTDBG("%s: Done splitting input %s, will clean up\n", __func__, input);
  fclose(inputFile);
  // need to close all output files
  for (int i = 0; i < numOutputFiles; i++)
  {
    fclose(outputFiles[i]);
  }
  free(outputFiles);
}

int main(int argc, char *argv[])
{
  // arg#1 -> input file name
  // arg#2 -> batch_size
  // args#3... -> output file names
  // flags: -b to use bytes (batch_size will be exact number of bytes instead of approximating to the closest line)
  if (argc < 4)
  {
    // TODO: document -r flag
    fprintf(stderr,
            "\n"
            "Usage: %s [-b] [-r] input_file batch_size output_file_1 output_file_2 [output_file_3 ...]\n\n"
            "    -b: use bytes (batch_size will be exact number of bytes instead of approximating to the closest line)\n\n",
            argv[0]);
    exit(1);
  }
  bool useBytes = 0, offset = 0, raw = 0;
  size_t batchSize = 0;
  char **outputFileNames = NULL;
  char *inputFileName = NULL;
  for (int i = 1; i < argc; i++)
  {
    //check flags
    if (strcmp(argv[i], "-b") == 0)
    {
      useBytes = 1;
      offset += 1;
      continue;
    }

    if (strcmp(argv[i], "-r") == 0)
    {
      raw = 1;
      offset += 1;
      continue;
    }

    int truePos = i - offset; //tracks true position without any added flags
    if (truePos == 1)
    {
      inputFileName = calloc(strlen(argv[i]) + 1, sizeof(char));
      strcpy(inputFileName, argv[i]);
    }

    if (truePos == 2)
      batchSize = atoi(argv[i]); //if 0 r_split default to approximate value

    if (truePos == 3)
    {
      outputFileNames = argv + i;
    }
  }

  SplitInput(inputFileName, batchSize, outputFileNames, argc - offset - 3, useBytes, raw);

  PRINTDBG("SplitInput is done\n");
  return 0;
}
