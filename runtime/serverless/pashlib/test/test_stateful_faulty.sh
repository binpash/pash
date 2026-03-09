#!/bin/bash
cd $(dirname "$0")

JOB_ID="test-ft-sort"
FOLDER_ID="ft-scripts"

export LANG=C.UTF-8
export LC_ALL=C.UTF-8

cleanup() {
    trap - EXIT INT TERM
    kill -- -$$ 2>/dev/null || true
}
trap cleanup EXIT INT TERM

aws s3 cp scripts/sort-lambda-1-faulty.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/sort-lambda-1-faulty.sh
aws s3 cp scripts/sort-lambda-2-faulty.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/sort-lambda-2-faulty.sh

SCRIPT_ID_1="sort-lambda-1-faulty"
SCRIPT_ID_2="sort-lambda-2-faulty"

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_1}\"],
    \"is_stateless\": false
    }" \
    --invocation-type Event \
    response.json

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_2}\"],
    \"is_stateless\": false
    }" \
    --invocation-type Event \
    response.json

bash scripts/sort-ec2-faulty.sh

faultless_out="/tmp/leash-ft-stateful-faultless.out"
rm -f ${faultless_out}
aws s3 cp s3://${AWS_BUCKET}/ft/stateful-faulty.txt  ${faultless_out}

echo "Finished running fault-free execution"

echo "================================"

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_1}\"],
    \"is_stateless\": false,
    \"time_to_crash\": 5
    }" \
    --invocation-type Event \
    response.json

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_2}\"],
    \"is_stateless\": false
    }" \
    --invocation-type Event \
    response.json

bash scripts/sort-ec2-faulty.sh

echo "Finished execution, checking results..."
faulty_out="/tmp/leash-ft-stateful-faulty.out"
rm -f ${faulty_out}
aws s3 cp s3://${AWS_BUCKET}/ft/stateful-faulty.txt  ${faulty_out}
echo "Comparing recovery vs fault-free execution..."
if diff -q ${faultless_out} ${faulty_out} > /dev/null; then
    echo "✅ Successfully recovered from failure, outputs match!"
else
    echo "❌ Recovery failed, outputs do NOT match!"
fi

rm -f ${faultless_out}
