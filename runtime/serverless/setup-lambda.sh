#!/usr/bin/env bash
set -euo pipefail

# Creates/updates the Lambda function using AWS CLI from the *current folder*.
# Expects your folder to contain: lambda-function.py, runtime/, aws/, etc.
#
# Required:
#   - AWS CLI configured (aws sts get-caller-identity works)
# Optional env vars:
#   - AWS_ACCOUNT_ID (if not set, we auto-detect)
#   - AWS_BUCKET (if not set, we leave it empty in Lambda env)
#   - ROLE_NAME (default: pash-lambda-exec)
#   - FUNCTION_NAME (default: lambda)
#   - REGION (default: us-east-1)

cd "$(dirname "$0")"

REGION="${REGION:-us-east-1}"
FUNCTION_NAME="${FUNCTION_NAME:-lambda-auto}"
ROLE_NAME="${ROLE_NAME:-leash-lambda-exec}"

RUNTIME="python3.12"
HANDLER="lambda-function.lambda_handler"
MEMORY_SIZE="1769"
TIMEOUT="30"

# Detect account ID if not provided
if [ -z "$AWS_BUCKET" ]; then
  echo "Error: Unable to detect AWS bucket name. Please set AWS_ACCOUNT_ID env var."
  exit 1
fi
if [ -z "$AWS_ACCOUNT_ID" ]; then
  echo "Error: Unable to detect AWS account ID. Please set AWS_ACCOUNT_ID env var."
  exit 1
fi

ROLE_ARN="arn:aws:iam::${AWS_ACCOUNT_ID}:role/${ROLE_NAME}"

ZIP_NAME="lambda.zip"

echo "[1/6] Region: ${REGION}"
echo "[1/6] Account: ${AWS_ACCOUNT_ID}"
echo "[1/6] Bucket: ${AWS_BUCKET}"
echo "[1/6] Function: ${FUNCTION_NAME}"
echo "[1/6] Role: ${ROLE_NAME}"
echo

# ---- IAM role (create if missing) ----
echo "[2/6] Ensuring IAM role exists..."
if aws iam get-role --role-name "${ROLE_NAME}" >/dev/null 2>&1; then
  echo "  - Role exists: ${ROLE_NAME}"
else
  cat > /tmp/lambda-trust-policy.json <<'JSON'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": { "Service": "lambda.amazonaws.com" },
      "Action": "sts:AssumeRole"
    }
  ]
}
JSON

  aws iam create-role \
    --role-name "${ROLE_NAME}" \
    --assume-role-policy-document file:///tmp/lambda-trust-policy.json \
    >/dev/null

  echo "  - Created role: ${ROLE_NAME}"
fi

# Attach basic logging policy (safe to re-run)
aws iam attach-role-policy \
  --role-name "${ROLE_NAME}" \
  --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole \
  >/dev/null 2>&1 || true

# Inline policy roughly matching your serverless.yml broad permissions
echo "[3/6] Putting inline policy on role (broad permissions)..."
cat > /tmp/pash-inline-policy.json <<'JSON'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "PashBroadAccess",
      "Effect": "Allow",
      "Action": [
        "lambda:*",
        "logs:*",
        "cloudwatch:*",
        "dynamodb:*",
        "s3:*",
        "sqs:*",
        "iam:GetPolicy",
        "iam:GetPolicyVersion",
        "iam:GetRole",
        "iam:GetRolePolicy",
        "iam:ListAttachedRolePolicies",
        "iam:ListRolePolicies",
        "iam:ListRoles",
        "iam:PassRole"
      ],
      "Resource": "*"
    }
  ]
}
JSON

aws iam put-role-policy \
  --role-name "${ROLE_NAME}" \
  --policy-name "pash-inline" \
  --policy-document file:///tmp/pash-inline-policy.json \
  >/dev/null

# IAM propagation can take a moment; a short sleep helps avoid create-function failures.
sleep 6

# ---- Package zip (include runtime/** aws/** lambda-function.py ; exclude resumability/**) ----
echo "[4/6] Packaging zip: ${ZIP_NAME}"
rm -f "${ZIP_NAME}"

# We build the zip deterministically from the current directory.
# Exclude resumability/** and common junk.
zip -r "${ZIP_NAME}" \
  "lambda-function.py" \
  "runtime" \
  "aws" \
  -x "resumability/*" "resumability/**" \
  -x "*.pyc" "__pycache__/*" "__pycache__/**" ".DS_Store" \
  >/dev/null

echo "  - Zip size: $(du -h "${ZIP_NAME}" | awk '{print $1}')"
echo

# ---- Create or Update Lambda ----
echo "[5/6] Creating or updating Lambda function..."
if aws lambda get-function --function-name "${FUNCTION_NAME}" --region "${REGION}" >/dev/null 2>&1; then
  echo "  - Function exists, updating code + config..."

  aws lambda update-function-code \
    --function-name "${FUNCTION_NAME}" \
    --region "${REGION}" \
    --zip-file "fileb://${ZIP_NAME}" \
    >/dev/null

  # wait a moment for the code update to propagate before updating config
  sleep 5

  aws lambda update-function-configuration \
    --function-name "${FUNCTION_NAME}" \
    --region "${REGION}" \
    --runtime "${RUNTIME}" \
    --handler "${HANDLER}" \
    --timeout "${TIMEOUT}" \
    --memory-size "${MEMORY_SIZE}" \
    --environment "Variables={AWS_ACCOUNT_ID=${AWS_ACCOUNT_ID},AWS_BUCKET=${AWS_BUCKET}}" \
    >/dev/null

else
  echo "  - Function not found, creating..."

  aws lambda create-function \
    --function-name "${FUNCTION_NAME}" \
    --region "${REGION}" \
    --runtime "${RUNTIME}" \
    --role "${ROLE_ARN}" \
    --handler "${HANDLER}" \
    --timeout "${TIMEOUT}" \
    --memory-size "${MEMORY_SIZE}" \
    --environment "Variables={AWS_ACCOUNT_ID=${AWS_ACCOUNT_ID},AWS_BUCKET=${AWS_BUCKET}}" \
    --zip-file "fileb://${ZIP_NAME}" \
    >/dev/null
fi

# ---- Set async retry attempts to 0 (matches maximumRetryAttempts: 0) ----
echo "[6/6] Setting async maximumRetryAttempts=0..."
aws lambda put-function-event-invoke-config \
  --function-name "${FUNCTION_NAME}" \
  --region "${REGION}" \
  --maximum-retry-attempts 0 \
  >/dev/null

echo
echo "Cleaning up temporary files..."
rm -f "${ZIP_NAME}" /tmp/lambda-trust-policy.json /tmp/pash-inline-policy.json
echo "✅ Done."
echo "Function ARN:" $(
    aws lambda get-function-configuration \
        --function-name "${FUNCTION_NAME}" \
        --region "${REGION}" \
        --query FunctionArn \
        --output text
)
