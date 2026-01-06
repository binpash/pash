#!/bin/bash
set -e

# Worker 2 for job job-20251218-221625
echo "[Worker 2] Starting..."

# Environment setup
export AWS_BUCKET="inout741448956691"
export PASHLIB_PATH="/opt/pashlib"

# Create FIFO for output (remove if exists from previous run)
FIFO="/tmp/split_2_output"
rm -f "$FIFO"
mkfifo "$FIFO"

# Background S3 upload using existing s3-put-object.py
# Usage: s3-put-object.py <object_key> <infile>
python3.9 aws/s3-put-object.py "jobs/job-20251218-221625/split_2.txt" "$FIFO" &
S3_PID=$!

# Run s3-shard-reader.py (at top level)
python3.9 aws/s3-shard-reader-direct.py \
    "oneliners/inputs/1M.txt" \
    "$FIFO" \
    "bytes=524288-786431" \
    split=2 \
    num_splits=4 \
    job_uid=job-20251218-221625

# Wait for S3 upload to complete
wait $S3_PID

# Cleanup
rm -f "$FIFO"

echo "[Worker 2] Complete: jobs/job-20251218-221625/split_2.txt"
