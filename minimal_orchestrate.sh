#!/bin/bash
# Simple orchestration for minimal ec2-lambda holepunch test
# This runs:
#   - EC2 script locally (sends data)
#   - Lambda script remotely (receives data)

set -e

export PASH_TOP="/home/ubuntu/pash"
INVOKE_SCRIPT="$PASH_TOP/aws/invoke-lambda.py"

# Use timestamp as folder ID for unique S3 paths
FOLDER_ID=$(date +%s)

# Script identifier for Lambda (can be any unique name)
LAMBDA_SCRIPT_ID="lambda"
LAMBDA2_SCRIPT_ID="lambda2"

echo "=========================================="
echo "Minimal EC2-Lambda Holepunch Test"
echo "=========================================="
echo "Folder ID: $FOLDER_ID"
echo "S3 Bucket: $AWS_BUCKET"
echo "Connection ID: 42"
echo "=========================================="
echo

echo "STEP 0: Delete RDV entry for connection ID '42' (if exists)..."
aws dynamodb delete-item \
  --table-name rdv \
  --key '{"key": {"S": "42"}}'
echo "  ✓ RDV entry deleted (if it existed)"
echo

aws dynamodb delete-item \
  --table-name rdv \
  --key '{"key": {"S": "43"}}'
echo "  ✓ RDV entry2 deleted (if it existed)"
echo

# Step 1: Upload Lambda script to S3
echo "[Step 1] Uploading Lambda script to S3..."
aws s3 cp "$PASH_TOP/minimal/lambda.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA_SCRIPT_ID}.sh"
echo "  ✓ Lambda script uploaded"
echo

echo "[Step 1b] Uploading Lambda2 script to S3..."
aws s3 cp "$PASH_TOP/minimal/lambda2.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA2_SCRIPT_ID}.sh"
echo "  ✓ Lambda2 script uploaded"
echo

# Step 2: Start both nodes in parallel (CRITICAL for holepunch timing!)
echo "[Step 2] Starting nodes in parallel..."
echo "  - EC2 (LOCAL): sending data via connection 42"

START_TIME=$(date +%s)

(cd "$PASH_TOP" && /bin/bash "$PASH_TOP/minimal/ec2.sh" "" "$FOLDER_ID") &
EC2_PID=$!

echo "  - Lambda (REMOTE): receiving data via connection 42"
python3 "$INVOKE_SCRIPT" "$LAMBDA_SCRIPT_ID" "$FOLDER_ID" "lambda" &
LAMBDA_PID=$!


echo "  - Lambda2 (REMOTE): receiving data via connection 43"
python3 "$INVOKE_SCRIPT" "$LAMBDA2_SCRIPT_ID" "$FOLDER_ID" "lambda" &
LAMBDA2_PID=$!

echo "  - Waiting for both nodes to complete..."
wait $EC2_PID $LAMBDA_PID $LAMBDA2_PID
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
    echo "✓ Holepunch test completed successfully!"
    echo ""
    echo "What happened:"
    echo "  1. EC2 got its STUN address and registered in DynamoDB"
    echo "  2. Lambda got its STUN address and registered in DynamoDB"
    echo "  3. Both found each other via connection ID '42'"
    echo "  4. TCP holepunch succeeded (simultaneous open)"
    echo "  5. Data transferred: EC2 → Lambda via holepunch"
else
    echo "✗ Holepunch test failed with exit code: $EXIT_CODE"
    echo ""
    echo "Debugging steps:"
    echo "  1. Check Lambda logs: aws logs tail /aws/lambda/lambda --follow"
    echo "  2. Check DynamoDB: aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"42\"}}'"
    echo "  3. Verify both nodes registered (should see N_0 and N_1 attributes)"
    echo "  4. Check timing - Lambda cold start may have taken too long"
fi
echo "=========================================="

exit $EXIT_CODE
