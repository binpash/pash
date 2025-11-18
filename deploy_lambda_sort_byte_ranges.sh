#!/bin/bash
# Deploy lambda-sort-ranges worker function (byte range version)

set -e

echo "=========================================="
echo "Deploying Lambda Sort Worker (Byte Ranges)"
echo "=========================================="

# Check environment variables
if [ -z "$AWS_ACCOUNT_ID" ]; then
    echo "Error: AWS_ACCOUNT_ID environment variable not set"
    exit 1
fi

# Package the Lambda function
echo "[Step 1] Packaging Lambda function..."
cd /home/ubuntu/pash
zip -q lambda_sort_worker_byte_ranges.zip lambda_sort_worker_byte_ranges.py
echo "  ✓ Created lambda_sort_worker_byte_ranges.zip"

# Check if function exists
FUNCTION_EXISTS=$(aws lambda get-function --function-name lambda-sort-ranges 2>&1 | grep -c "ResourceNotFoundException" || true)

if [ "$FUNCTION_EXISTS" -eq 1 ]; then
    # Create new function
    echo "[Step 2] Creating new Lambda function 'lambda-sort-ranges'..."

    aws lambda create-function \
        --function-name lambda-sort-ranges \
        --runtime python3.9 \
        --role arn:aws:iam::${AWS_ACCOUNT_ID}:role/pash-release-us-east-1-lambdaRole \
        --handler lambda_sort_worker_byte_ranges.lambda_handler \
        --zip-file fileb://lambda_sort_worker_byte_ranges.zip \
        --timeout 300 \
        --memory-size 3008 \
        --ephemeral-storage Size=2048 \
        --region us-east-1

    echo "  ✓ Lambda function created"
else
    # Update existing function
    echo "[Step 2] Updating existing Lambda function 'lambda-sort-ranges'..."

    aws lambda update-function-code \
        --function-name lambda-sort-ranges \
        --zip-file fileb://lambda_sort_worker_byte_ranges.zip \
        --region us-east-1

    echo "  ✓ Lambda function updated"
fi

echo
echo "=========================================="
echo "✓ Deployment complete!"
echo "Function name: lambda-sort-ranges"
echo "=========================================="
echo
echo "You can now test with:"
echo "  python3 manual_s3_orchestrator_byte_ranges.py \\"
echo "    --bucket \$AWS_BUCKET \\"
echo "    --input oneliners/inputs/1M.txt \\"
echo "    --output oneliners/outputs/byte-range-result.txt \\"
echo "    --workers 2"
