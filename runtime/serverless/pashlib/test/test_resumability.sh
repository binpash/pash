#!/bin/bash
cd $(dirname "$0")

JOB_ID="test-resumability"
FOLDER_ID="resumability-scripts"

export LANG=C.UTF-8
export LC_ALL=C.UTF-8

cleanup() {
    trap - EXIT INT TERM
    kill -- -$$ 2>/dev/null || true
}
trap cleanup EXIT INT TERM

aws s3 cp scripts/resumability-lambda.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/resumability-lambda.sh
SCRIPT_ID_1="resumability-lambda"

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_1}\"],
    \"is_stateless\": false,
    \"timeout\": 900
    }" \
    --invocation-type Event \
    response.json

bash scripts/resumability-ec2.sh

no_resumability_trigger="/tmp/leash-ft-no-resumability-triggered.out"
rm -f ${no_resumability_trigger}
aws s3 cp s3://${AWS_BUCKET}/resumability/resumability-sorted.txt  ${no_resumability_trigger}

echo "Finished baseline execution"

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_1}\"],
    \"is_stateless\": false,
    \"timeout\": 30
    }" \
    --invocation-type Event \
    response.json

bash scripts/resumability-ec2.sh

resumability_trigger="/tmp/leash-ft-resumability-triggered.out"
rm -f ${resumability_trigger}
aws s3 cp s3://${AWS_BUCKET}/resumability/resumability-sorted.txt  ${resumability_trigger}  

echo "Finished resumability execution, checking results..."

if diff -q ${no_resumability_trigger} ${resumability_trigger} > /dev/null; then
    echo "✅ Successfully resumed from timeout, outputs match!"
else
    echo "❌ Resumability failed, outputs do NOT match!"
fi