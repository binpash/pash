#!/bin/bash
# Hybrid Lambda-EC2 orchestration for test_merge2cat
# Lambda workers: l1.sh, l2.sh, sortl3.sh, sortl4.sh
# EC2 workers: merger1.sh, sortmerge.sh

set -e

export PASH_TOP="/home/ubuntu/pash"
INVOKE_SCRIPT="$PASH_TOP/aws/invoke-lambda.py"

# Timestamp-based folder ID for unique S3 paths
FOLDER_ID=$(date +%s)

# Script identifiers for Lambda
L1_SCRIPT_ID="l1"
L2_SCRIPT_ID="l2"
SORTL3_SCRIPT_ID="sortl3"
SORTL4_SCRIPT_ID="sortl4"

echo "=========================================="
echo "Hybrid Lambda-EC2 Test (Merge2Cat)"
echo "=========================================="
echo "Folder ID: $FOLDER_ID"
echo "S3 Bucket: $AWS_BUCKET"
echo "Lambda: l1, l2, sortl3, sortl4"
echo "EC2: merger1, sortmerge"
echo "Connections: aa, bb, cc, dd, ee, ff"
echo "=========================================="
echo

# STEP 0: Clean RDV entries for all 6 connections
echo "STEP 0: Cleaning RDV entries..."

aws dynamodb delete-item --table-name rdv --key '{"key": {"S": "aa"}}' 2>/dev/null || true
echo "  ✓ Deleted RDV entry: aa (L1→Merger1)"

aws dynamodb delete-item --table-name rdv --key '{"key": {"S": "bb"}}' 2>/dev/null || true
echo "  ✓ Deleted RDV entry: bb (L2→Merger1)"

aws dynamodb delete-item --table-name rdv --key '{"key": {"S": "cc"}}' 2>/dev/null || true
echo "  ✓ Deleted RDV entry: cc (Merger1→Sortl3)"

aws dynamodb delete-item --table-name rdv --key '{"key": {"S": "dd"}}' 2>/dev/null || true
echo "  ✓ Deleted RDV entry: dd (Sortl3→Sortmerge)"

aws dynamodb delete-item --table-name rdv --key '{"key": {"S": "ee"}}' 2>/dev/null || true
echo "  ✓ Deleted RDV entry: ee (Merger1→Sortl4)"

aws dynamodb delete-item --table-name rdv --key '{"key": {"S": "ff"}}' 2>/dev/null || true
echo "  ✓ Deleted RDV entry: ff (Sortl4→Sortmerge)"

# Clean job_uid-based RDV entries (for S3 shard reader inter-shard communication)
aws dynamodb delete-item --table-name rdv --key '{"key": {"S": "j1-split-0-to-1"}}' 2>/dev/null || true
echo "  ✓ Deleted RDV entry: j1 (job1 shard coordination)"

echo "  ✓ All RDV entries cleaned"
echo

# STEP 1: Upload Lambda scripts to S3
echo "[Step 1] Uploading Lambda scripts to S3..."

aws s3 cp "$PASH_TOP/testing/test_merge2cat/l1.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${L1_SCRIPT_ID}.sh"
echo "  ✓ Uploaded l1.sh"

aws s3 cp "$PASH_TOP/testing/test_merge2cat/l2.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${L2_SCRIPT_ID}.sh"
echo "  ✓ Uploaded l2.sh"

aws s3 cp "$PASH_TOP/testing/test_merge2cat/sortl3.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${SORTL3_SCRIPT_ID}.sh"
echo "  ✓ Uploaded sortl3.sh"

aws s3 cp "$PASH_TOP/testing/test_merge2cat/sortl4.sh" \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${SORTL4_SCRIPT_ID}.sh"
echo "  ✓ Uploaded sortl4.sh"

echo "  ✓ All Lambda scripts uploaded"
echo

# STEP 2: Start all 6 components in parallel (CRITICAL for holepunch!)
echo "[Step 2] Starting all components in parallel..."
START_TIME=$(date +%s)

# Lambda workers (4 total)
echo "  - L1 (LAMBDA): reading S3 shard 0, sending to Merger1 (conn: aa)"
python3 "$INVOKE_SCRIPT" "$L1_SCRIPT_ID" "$FOLDER_ID" "lambda" &
L1_PID=$!

echo "  - L2 (LAMBDA): reading S3 shard 1, sending to Merger1 (conn: bb)"
python3 "$INVOKE_SCRIPT" "$L2_SCRIPT_ID" "$FOLDER_ID" "lambda" &
L2_PID=$!

echo "  - Sortl3 (LAMBDA): receiving from Merger1 (conn: cc), sorting, sending to Sortmerge (conn: dd)"
python3 "$INVOKE_SCRIPT" "$SORTL3_SCRIPT_ID" "$FOLDER_ID" "lambda" &
SORTL3_PID=$!

echo "  - Sortl4 (LAMBDA): receiving from Merger1 (conn: ee), sorting, sending to Sortmerge (conn: ff)"
python3 "$INVOKE_SCRIPT" "$SORTL4_SCRIPT_ID" "$FOLDER_ID" "lambda" &
SORTL4_PID=$!

# EC2 workers (2 total)
echo "  - Merger1 (EC2): receiving from L1+L2 (conn: aa+bb), merging, splitting to Sortl3+Sortl4 (conn: cc+ee)"
(cd "$PASH_TOP" && /bin/bash testing/test_merge2cat/merger1.sh "" "$FOLDER_ID") &
MERGER1_PID=$!

echo "  - Sortmerge (EC2): receiving from Sortl3+Sortl4 (conn: dd+ff), merging, writing to S3"
(cd "$PASH_TOP" && /bin/bash testing/test_merge2cat/sortmerge.sh "" "$FOLDER_ID") &
SORTMERGE_PID=$!

echo ""
echo "All 6 components launched. Waiting for completion..."
wait $L1_PID $L2_PID $SORTL3_PID $SORTL4_PID $MERGER1_PID $SORTMERGE_PID
EXIT_CODE=$?

END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Hybrid Lambda-EC2 test completed successfully!"
    echo ""
    echo "Data flow:"
    echo "  1. L1 (Lambda) read S3 shard 0 → Merger1 (EC2) via aa"
    echo "  2. L2 (Lambda) read S3 shard 1 → Merger1 (EC2) via bb"
    echo "  3. Merger1 (EC2) merged+split → Sortl3 (Lambda) via cc + Sortl4 (Lambda) via ee"
    echo "  4. Sortl3 (Lambda) sorted → Sortmerge (EC2) via dd"
    echo "  5. Sortl4 (Lambda) sorted → Sortmerge (EC2) via ff"
    echo "  6. Sortmerge (EC2) final merge → S3 output"
    echo ""
    echo "Total execution time: ${ELAPSED_TIME}s"
else
    echo "✗ Hybrid test failed with exit code: $EXIT_CODE"
    echo ""
    echo "Debugging steps:"
    echo "  1. Check Lambda logs:"
    echo "     aws logs tail /aws/lambda/lambda --follow"
    echo "  2. Check RDV entries for connection state:"
    echo "     aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"aa\"}}'"
    echo "     aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"bb\"}}'"
    echo "     aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"cc\"}}'"
    echo "     aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"dd\"}}'"
    echo "     aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"ee\"}}'"
    echo "     aws dynamodb get-item --table-name rdv --key '{\"key\":{\"S\":\"ff\"}}'"
    echo "  3. Verify all nodes registered (should see N_0 and N_1 attributes)"
    echo "  4. Check Lambda cold start timing (may cause holepunch timeout)"
fi
echo "=========================================="

exit $EXIT_CODE
