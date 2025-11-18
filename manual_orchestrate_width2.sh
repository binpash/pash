#!/bin/bash
# Manual orchestration for width=2 sort using pre-generated scripts
# This mimics how pa.sh + ServerlessManager execute scripts:
#   - EC2 scripts (split/merge) run LOCALLY
#   - Lambda scripts run on AWS Lambda

set -e

SCRIPTS_DIR="/home/ubuntu/pash/evaluation/benchmarks/oneliners/logs/sort.sh:1M.txt:2/scripts/1762478811"
export PASH_TOP="/home/ubuntu/pash"
INVOKE_SCRIPT="$PASH_TOP/aws/invoke-lambda.py"

# Script IDs from pre-generated logs
SPLITTER="559c679d-49ed-49c3-9a1f-08c3ce9e5f02"  # EC2: download + split
LAMBDA1="7422cb0d-628b-49ae-b29d-387640f949f6"   # Lambda: sort
LAMBDA2="6b9eae70-3757-42a4-a3fb-eee1dc6a73ca"   # Lambda: sort
MERGER="839ba2b4-ed80-4864-a800-04694e745f76"    # EC2: merge + upload

# Use timestamp as folder ID
FOLDER_ID=$(date +%s)

echo "=========================================="
echo "Manual Width=2 Sort Orchestration"
echo "=========================================="
echo "Folder ID: $FOLDER_ID"
echo "Scripts: $SCRIPTS_DIR"
echo "S3 Bucket: $AWS_BUCKET"
echo "=========================================="
echo

# Step 1: Prepare scripts
echo "[Step 1] Preparing scripts..."

# Copy EC2 scripts to /tmp (run locally like pa.sh does)
echo "  - Copying EC2 scripts to /tmp for local execution..."
cp "$SCRIPTS_DIR/${SPLITTER}.sh" "/tmp/$SPLITTER"
cp "$SCRIPTS_DIR/${MERGER}.sh" "/tmp/$MERGER"
chmod +x "/tmp/$SPLITTER" "/tmp/$MERGER"

# Upload Lambda scripts to S3
echo "  - Uploading Lambda scripts to S3..."
aws s3 cp "$SCRIPTS_DIR/${LAMBDA1}.sh" "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA1}.sh"
aws s3 cp "$SCRIPTS_DIR/${LAMBDA2}.sh" "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA2}.sh"
echo "  ✓ Scripts prepared"
echo

# Step 2: Execute all nodes in parallel (critical for NAT traversal timing)
# EC2 scripts run LOCALLY, Lambda scripts run on AWS Lambda
echo "[Step 2] Starting all nodes in parallel..."
echo "  - Splitter (LOCAL EC2)"
(cd "$PASH_TOP" && /bin/bash "/tmp/$SPLITTER" "$FOLDER_ID") &
SPLITTER_PID=$!

echo "  - Lambda Worker 1"
python3 "$INVOKE_SCRIPT" "$LAMBDA1" "$FOLDER_ID" "lambda" &
LAMBDA1_PID=$!

echo "  - Lambda Worker 2"
python3 "$INVOKE_SCRIPT" "$LAMBDA2" "$FOLDER_ID" "lambda" &
LAMBDA2_PID=$!

echo "  - Merger (LOCAL EC2)"
(cd "$PASH_TOP" && /bin/bash "/tmp/$MERGER" "$FOLDER_ID") &
MERGER_PID=$!

echo "  - Waiting for all nodes to complete..."
wait $SPLITTER_PID $LAMBDA1_PID $LAMBDA2_PID $MERGER_PID
EXIT_CODE=$?
echo "  ✓ All nodes completed (exit code: $EXIT_CODE)"
echo

echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Orchestration completed successfully!"
else
    echo "✗ Orchestration failed with exit code: $EXIT_CODE"
fi
echo "Output: s3://$AWS_BUCKET/oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt"
echo "=========================================="

# Cleanup
echo
echo "Cleaning up temporary files..."
rm -f "/tmp/$SPLITTER" "/tmp/$MERGER"
echo "✓ Cleanup complete"
