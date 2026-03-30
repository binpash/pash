#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
SERVERLESS_DIR="$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)"

REGION="${REGION:-${AWS_REGION:-${AWS_DEFAULT_REGION:-us-east-1}}}"
FUNCTION_NAME="${FUNCTION_NAME:-lambda}"
ROLE_NAME="${ROLE_NAME:-leash-lambda-exec}"
REPOSITORY_NAME="${REPOSITORY_NAME:-pash-serverless}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
LOCAL_IMAGE_TAG="${LOCAL_IMAGE_TAG:-${REPOSITORY_NAME}:${IMAGE_TAG}}"
TIMEOUT="${1:-900}"
MEMORY_SIZE="${MEMORY_SIZE:-1769}"
EPHEMERAL_DISK_SIZE="${EPHEMERAL_DISK_SIZE:-10240}"
DOCKER_BIN="${DOCKER_BIN:-sudo docker}"

cd "$SERVERLESS_DIR"

require_env() {
  local var_name="$1"
  if [[ -z "${!var_name:-}" ]]; then
    echo "Error: Please export ${var_name} in your shell before running this script." >&2
    exit 1
  fi
}

require_env AWS_ACCOUNT_ID
require_env AWS_BUCKET

ROLE_ARN="arn:aws:iam::${AWS_ACCOUNT_ID}:role/${ROLE_NAME}"
ECR_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${REPOSITORY_NAME}:${IMAGE_TAG}"
BASIC_EXEC_POLICY_ARN="arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole"
INLINE_POLICY_NAME="pash-inline"

echo "[1/7] Region: ${REGION}"
echo "[1/7] Account: ${AWS_ACCOUNT_ID}"
echo "[1/7] Bucket: ${AWS_BUCKET}"
echo "[1/7] Function: ${FUNCTION_NAME}"
echo "[1/7] Role: ${ROLE_NAME}"
echo "[1/7] Repository: ${REPOSITORY_NAME}"
echo "[1/7] Local image: ${LOCAL_IMAGE_TAG}"
echo "[1/7] Remote image: ${ECR_URI}"
echo

echo "[2/7] Ensuring IAM role exists..."
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

policy_updated=0
attached_policy_count="$(
  aws iam list-attached-role-policies \
    --role-name "${ROLE_NAME}" \
    --query "length(AttachedPolicies[?PolicyArn=='${BASIC_EXEC_POLICY_ARN}'])" \
    --output text 2>/dev/null || echo "0"
)"

if [[ "${attached_policy_count}" == "0" ]]; then
  aws iam attach-role-policy \
    --role-name "${ROLE_NAME}" \
    --policy-arn "${BASIC_EXEC_POLICY_ARN}" \
    >/dev/null
  policy_updated=1
  echo "  - Attached managed policy: AWSLambdaBasicExecutionRole"
else
  echo "  - Managed policy already attached: AWSLambdaBasicExecutionRole"
fi

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
        "ecr:GetAuthorizationToken",
        "ecr:BatchCheckLayerAvailability",
        "ecr:CompleteLayerUpload",
        "ecr:CreateRepository",
        "ecr:DescribeRepositories",
        "ecr:InitiateLayerUpload",
        "ecr:PutImage",
        "ecr:UploadLayerPart",
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

if aws iam get-role-policy \
  --role-name "${ROLE_NAME}" \
  --policy-name "${INLINE_POLICY_NAME}" \
  >/dev/null 2>&1; then
  echo "  - Inline policy already attached: ${INLINE_POLICY_NAME}"
else
  aws iam put-role-policy \
    --role-name "${ROLE_NAME}" \
    --policy-name "${INLINE_POLICY_NAME}" \
    --policy-document file:///tmp/pash-inline-policy.json \
    >/dev/null
  policy_updated=1
  echo "  - Attached inline policy: ${INLINE_POLICY_NAME}"
fi

if [[ "${policy_updated}" -eq 1 ]]; then
  sleep 6
else
  echo "  - Skipping IAM propagation sleep (no role policy changes)."
fi

echo "[3/7] Ensuring ECR repository exists..."
if aws ecr describe-repositories --repository-names "${REPOSITORY_NAME}" --region "${REGION}" >/dev/null 2>&1; then
  echo "  - Repository exists: ${REPOSITORY_NAME}"
else
  aws ecr create-repository --repository-name "${REPOSITORY_NAME}" --region "${REGION}" >/dev/null
  echo "  - Created repository: ${REPOSITORY_NAME}"
fi

echo "[4/7] Logging Docker into ECR..."
aws ecr get-login-password --region "${REGION}" | ${DOCKER_BIN} login --username AWS --password-stdin "${AWS_ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"

echo "[5/7] Tagging and pushing image..."
${DOCKER_BIN} tag "${LOCAL_IMAGE_TAG}" "${ECR_URI}"
${DOCKER_BIN} push "${ECR_URI}"

echo "[6/7] Creating or updating Lambda function..."
if aws lambda get-function --function-name "${FUNCTION_NAME}" --region "${REGION}" >/dev/null 2>&1; then
  echo "  - Function exists, updating image + config..."
  aws lambda update-function-code \
    --function-name "${FUNCTION_NAME}" \
    --region "${REGION}" \
    --image-uri "${ECR_URI}" \
    >/dev/null

  echo "  - Waiting for image update to finish before applying configuration..."
  aws lambda wait function-updated \
    --function-name "${FUNCTION_NAME}" \
    --region "${REGION}"

  aws lambda update-function-configuration \
    --function-name "${FUNCTION_NAME}" \
    --region "${REGION}" \
    --timeout "${TIMEOUT}" \
    --memory-size "${MEMORY_SIZE}" \
    --ephemeral-storage "Size=${EPHEMERAL_DISK_SIZE}" \
    --role "${ROLE_ARN}" \
    --environment "Variables={AWS_ACCOUNT_ID=${AWS_ACCOUNT_ID},AWS_BUCKET=${AWS_BUCKET}}" \
    >/dev/null
else
  echo "  - Function not found, creating..."
  aws lambda create-function \
    --function-name "${FUNCTION_NAME}" \
    --region "${REGION}" \
    --package-type Image \
    --code "ImageUri=${ECR_URI}" \
    --role "${ROLE_ARN}" \
    --timeout "${TIMEOUT}" \
    --memory-size "${MEMORY_SIZE}" \
    --ephemeral-storage "Size=${EPHEMERAL_DISK_SIZE}" \
    --environment "Variables={AWS_ACCOUNT_ID=${AWS_ACCOUNT_ID},AWS_BUCKET=${AWS_BUCKET}}" \
    >/dev/null
fi

echo "[7/7] Setting async maximumRetryAttempts=0..."
aws lambda put-function-event-invoke-config \
  --function-name "${FUNCTION_NAME}" \
  --region "${REGION}" \
  --maximum-retry-attempts 0 \
  >/dev/null

echo
echo "Cleaning up temporary files..."
rm -f /tmp/lambda-trust-policy.json /tmp/pash-inline-policy.json
