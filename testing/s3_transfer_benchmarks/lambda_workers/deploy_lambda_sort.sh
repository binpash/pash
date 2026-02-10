#!/usr/bin/env bash
# Deploy lambda-sort worker function

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Deploying Lambda Sort Worker"
echo "=========================================="

# Check environment variables
if [ -z "$AWS_ACCOUNT_ID" ]; then
    echo "Error: AWS_ACCOUNT_ID environment variable not set"
    exit 1
fi

# Package the Lambda function
echo "[Step 1] Packaging Lambda function..."
cd "$SCRIPT_DIR"
zip -q lambda_sort_worker.zip lambda_sort_worker.py
echo "  ✓ Created lambda_sort_worker.zip"

# Check if function exists
FUNCTION_EXISTS=$(aws lambda get-function --function-name lambda-sort 2>&1 | grep -c "ResourceNotFoundException" || true)

if [ "$FUNCTION_EXISTS" -eq 1 ]; then
    # Create new function
    echo "[Step 2] Creating new Lambda function 'lambda-sort'..."

    aws lambda create-function \
        --function-name lambda-sort \
        --runtime python3.9 \
        --role arn:aws:iam::${AWS_ACCOUNT_ID}:role/pash-release-us-east-1-lambdaRole \
        --handler lambda_sort_worker.lambda_handler \
        --zip-file fileb://lambda_sort_worker.zip \
        --timeout 300 \
        --memory-size 3008 \
        --ephemeral-storage Size=2048 \
        --region us-east-1

    echo "  ✓ Lambda function created"
else
    # Update existing function
    echo "[Step 2] Updating existing Lambda function 'lambda-sort'..."

    aws lambda update-function-code \
        --function-name lambda-sort \
        --zip-file fileb://lambda_sort_worker.zip \
        --region us-east-1

    echo "  ✓ Lambda function updated"
fi

echo
echo "=========================================="
echo "✓ Deployment complete!"
echo "Function name: lambda-sort"
echo "=========================================="
echo
echo "You can now test with:"
echo "  python3 ../orchestrators/manual_s3_orchestrator.py \\"
echo "    --bucket \$AWS_BUCKET \\"
echo "    --input oneliners/inputs/1M.txt \\"
echo "    --output oneliners/outputs/manual-sort-result.txt \\"
echo "    --workers 2"
