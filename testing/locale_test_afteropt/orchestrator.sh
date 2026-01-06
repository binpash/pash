#!/bin/bash
# Distributed sort orchestration test
# This runs:
#   - splitter.sh locally (reads S3, splits data to connections 34 and 2b)
#   - l1.sh on Lambda (receives from 34, sorts, sends to 17)
#   - l2.sh on Lambda (receives from 2b, sorts, sends to 4f)
#   - merger.sh locally (receives from 17 and 4f, merges, writes to S3)

set -e



export PASH_TOP="/home/ubuntu/pash"
INVOKE_SCRIPT="$PASH_TOP/aws/invoke-lambda.py"

# Use timestamp as folder ID for unique S3 paths
FOLDER_ID=$(date +%s)

# Script identifier for Lambda (can be any unique name)
LAMBDA_SCRIPT_ID="l1"
LAMBDA2_SCRIPT_ID="l2"
MERGER_SCRIPT_ID="merger"

echo "=========================================="
echo "Distributed Sort Test (Splitter-Lambda-Merger)"
echo "=========================================="
echo "Folder ID: $FOLDER_ID"
echo "S3 Bucket: $AWS_BUCKET"
echo "Connections: 34, 2b, 17, 4f"
echo "=========================================="
echo

# 34 17 2b 4f are the identifiers for the holepunch connections
# used in this test (hardcoded in the scripts)


echo "STEP 0: Delete RDV entry for connection ID '34' (if exists)..."
aws dynamodb delete-item \
  --table-name rdv \
  --key '{"key": {"S": "50"}}'
echo "  ✓ RDV entry deleted (if it existed)"
echo

echo "STEP 0: Delete RDV entry for connection ID '17' (if exists)..."
aws dynamodb delete-item \
  --table-name rdv \
  --key '{"key": {"S": "ed"}}'
echo "  ✓ RDV entry deleted (if it existed)"
echo

echo "STEP 0: Delete RDV entry for connection ID '19' (if exists)..."
aws dynamodb delete-item \
  --table-name rdv \
  --key '{"key": {"S": "19"}}'
echo "  ✓ RDV entry deleted (if it existed)"
echo



# Step 1: Upload Lambda script to S3
echo "[Step 1] Uploading Lambda script to S3..."
aws s3 cp "l1.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA_SCRIPT_ID}.sh"
echo "  ✓ Lambda script uploaded"
echo

echo "[Step 1b] Uploading Lambda2 script to S3..."
aws s3 cp "l2.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA2_SCRIPT_ID}.sh"
echo "  ✓ Lambda2 script uploaded"
echo

# Step 2: Start all 4 nodes in parallel (CRITICAL for holepunch timing!)
echo "[Step 2] Starting nodes in parallel..."

START_TIME=$(date +%s)

echo "  - Merger (LOCAL): receiving from connections ed and 4f, merging output"
(cd "$PASH_TOP" && /bin/bash testing/locale_test_afteropt/merger.sh "" "$FOLDER_ID") &
MERGER_PID=$!

echo "  - Lambda1 (REMOTE): receiving from connection 50, sorting, sending to ed"
python3 "$INVOKE_SCRIPT" "$LAMBDA_SCRIPT_ID" "$FOLDER_ID" "lambda" &
LAMBDA_PID=$!

echo "  - Lambda2 (REMOTE): receiving from connection ed, sorting, sending to 4f"
python3 "$INVOKE_SCRIPT" "$LAMBDA2_SCRIPT_ID" "$FOLDER_ID" "lambda" &
LAMBDA2_PID=$!


echo "  - Waiting for all 3 nodes to complete..."
wait $LAMBDA_PID $LAMBDA2_PID $MERGER_PID
EXIT_CODE=$?
echo "  ✓ All nodes completed (exit code: $EXIT_CODE)"
echo

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
echo "Total time taken for holepunch test: $ELAPSED_TIME seconds"
echo

# Results
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Distributed sort test completed successfully!"
    echo ""
    echo "What happened:"
    echo "  1. Splitter read from S3 and split data to connections 34 and 2b"
    echo "  2. Lambda1 received from 34, sorted, sent to 17"
    echo "  3. Lambda2 received from 2b, sorted, sent to 4f"
    echo "  4. Merger received from 17 and 4f, merged with sort -m"
    echo "  5. Final sorted output written to S3"
else
    echo "✗ Distributed sort test failed with exit code: $EXIT_CODE"
    echo ""
    echo "Debugging steps:"
    echo "  1. Check Lambda logs: aws logs tail /aws/lambda/lambda --follow"
    echo "  2. Check DynamoDB for connection 34: aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"34\"}}'"
    echo "  3. Check DynamoDB for connection 2b: aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"2b\"}}'"
    echo "  4. Check DynamoDB for connection 17: aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"17\"}}'"
    echo "  5. Check DynamoDB for connection 4f: aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"4f\"}}'"
    echo "  6. Verify all nodes registered (should see N_0 and N_1 attributes)"
    echo "  7. Check timing - Lambda cold start may have taken too long"
fi
echo "=========================================="

exit $EXIT_CODE
