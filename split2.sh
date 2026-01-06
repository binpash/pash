#!/usr/bin/env bash
set -euo pipefail

########################################
# Configuration
########################################

WORKDIR="$PWD/pipeline_io"
DEBUGDIR="$WORKDIR/debug"
PIPEDIR="$WORKDIR/pipes"


DEFAULT_S3_KEY="oneliners/inputs/1M.txt"
S3_KEY="${1:-$DEFAULT_S3_KEY}"


SPLIT_SIZE_BYTES=1000000

########################################
# Directory layout
########################################

mkdir -p "$DEBUGDIR" "$PIPEDIR"

########################################
# Named pipes (descriptive)
########################################

S3_STREAM_PIPE="$PIPEDIR/s3_object_stream.pipe"
SPLIT_INPUT_PIPE="$PIPEDIR/split_input.pipe"
SPLIT_FIRST_PIPE="$PIPEDIR/split_first_chunk.pipe"
SPLIT_REST_PIPE="$PIPEDIR/split_remaining.pipe"

########################################
# Cleanup helpers
########################################

cleanup() {
  rm -f \
    "$S3_STREAM_PIPE" \
    "$SPLIT_INPUT_PIPE" \
    "$SPLIT_FIRST_PIPE" \
    "$SPLIT_REST_PIPE"
}
trap cleanup EXIT INT TERM

########################################
# Recreate pipes
########################################

cleanup

mkfifo \
  "$S3_STREAM_PIPE" \
  "$SPLIT_INPUT_PIPE" \
  "$SPLIT_FIRST_PIPE" \
  "$SPLIT_REST_PIPE"

########################################
# Streaming pipeline
########################################

# 1) Stream object from S3
python3.9 aws/s3-get-object.py \
  "$S3_KEY" \
  "$S3_STREAM_PIPE" &
PID_S3=$!

# Save raw S3 stream
tee "$DEBUGDIR/s3_raw_stream.txt" < "$S3_STREAM_PIPE" > "$SPLIT_INPUT_PIPE" &
PID_TEE_S3=$!

# 2) Split stream
runtime/r_split -r \
  "$SPLIT_INPUT_PIPE" \
  "$SPLIT_SIZE_BYTES" \
  "$SPLIT_FIRST_PIPE" \
  "$SPLIT_REST_PIPE" &
PID_SPLIT=$!

# 3) Save split outputs
cat "$SPLIT_FIRST_PIPE" > "$DEBUGDIR/split_first_chunk.txt" &
PID_SAVE_FIRST=$!

cat "$SPLIT_REST_PIPE" > "$DEBUGDIR/split_remaining.txt" &
PID_SAVE_REST=$!

########################################
# Wait for completion
########################################

wait \
  "$PID_S3" \
  "$PID_TEE_S3" \
  "$PID_SPLIT" \
  "$PID_SAVE_FIRST" \
  "$PID_SAVE_REST"
