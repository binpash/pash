#!/bin/bash
# Merged local sort (no holepunch, all on EC2)
# Hardcoded byte ranges: shard0=0-524332, shard1=524333-1048575

set -e

export PASH_TOP="/home/ubuntu/pash"
export AWS_BUCKET="${AWS_BUCKET:-inout741448956691}"

FOLDER_ID=${2:-$(date +%s)}

# Create temp directory for fifos
TMPDIR="/tmp/pash_local_$$"
mkdir -p "$TMPDIR"
cd "$TMPDIR"

# Define FIFOs
FIFO_RAW0="raw0.fifo"
FIFO_RAW1="raw1.fifo"
FIFO_SORTED0="sorted0.fifo"
FIFO_SORTED1="sorted1.fifo"
FIFO_MERGED="merged.fifo"

# Cleanup function
cleanup() {
    rm -f "$FIFO_RAW0" "$FIFO_RAW1" "$FIFO_SORTED0" "$FIFO_SORTED1" "$FIFO_MERGED"
}

trap cleanup EXIT

# Create fifos
mkfifo "$FIFO_RAW0" "$FIFO_RAW1" "$FIFO_SORTED0" "$FIFO_SORTED1" "$FIFO_MERGED"

echo "=========================================="
echo "Local Parallel Sort (No Holepunch)"
echo "=========================================="
echo "S3 Input: s3://$AWS_BUCKET/oneliners/inputs/1M.txt"
echo "S3 Output: s3://$AWS_BUCKET/oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt"
echo "Folder ID: $FOLDER_ID"
echo

# Downloader 1: Download shard 0 (bytes 0-524332)
echo "[Downloader 1] Downloading bytes 0-524332..."
{
    python3.9 "$PASH_TOP/aws/s3-get-object.py" \
        "oneliners/inputs/1M.txt" \
        "$FIFO_RAW0" \
        "bytes=0-524332"
} &
PID_DL1=$!

# Downloader 2: Download shard 1 (bytes 524333-1048575)
echo "[Downloader 2] Downloading bytes 524333-1048575..."
{
    python3.9 "$PASH_TOP/aws/s3-get-object.py" \
        "oneliners/inputs/1M.txt" \
        "$FIFO_RAW1" \
        "bytes=524333-1048575"
} &
PID_DL2=$!

# Sorter 1: Sort shard 0
{
    sort < "$FIFO_RAW0" > "$FIFO_SORTED0"
    echo "[Sorter 1] Done sorting shard 0"
} &
PID_SORT1=$!

# Sorter 2: Sort shard 1
{
    sort < "$FIFO_RAW1" > "$FIFO_SORTED1"
    echo "[Sorter 2] Done sorting shard 1"
} &
PID_SORT2=$!

# Merger: Merge both sorted streams
echo "[Merger] Waiting for sorted shards to merge..."
{
    sort -m "$FIFO_SORTED0" "$FIFO_SORTED1" > "$FIFO_MERGED"
} &
PID_MERGE=$!

# Uploader: Upload merged result to S3
echo "[Uploader] Waiting for merged data..."
{
    python3.9 "$PASH_TOP/aws/s3-put-object.py" \
        "oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt" \
        "$FIFO_MERGED"
    echo "[Uploader] Upload complete"
} &
PID_UPLOAD=$!

# Wait for all background processes
echo
echo "Waiting for pipeline to complete..."
wait $PID_DL1 $PID_DL2 $PID_SORT1 $PID_SORT2 $PID_MERGE $PID_UPLOAD

EXIT_CODE=$?

echo
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Local parallel sort completed successfully!"
else
    echo "✗ Pipeline failed with exit code: $EXIT_CODE"
fi
echo "=========================================="
echo

exit $EXIT_CODE
