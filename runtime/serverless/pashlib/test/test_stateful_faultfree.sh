#!/bin/bash
cd $(dirname "$0")

JOB_ID="test-ft-sort"
FOLDER_ID="ft-scripts"

export LANG=C.UTF-8
export LC_ALL=C.UTF-8

aws s3 cp scripts/sort-lambda-1.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/sort-lambda-1.sh
aws s3 cp scripts/sort-lambda-2.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/sort-lambda-2.sh
SCRIPT_ID_1="sort-lambda-1"
SCRIPT_ID_2="sort-lambda-2"

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

bash scripts/sort-ec2.sh

echo "Finished execution, checking results..."

input_file="/tmp/1M.txt"
base_line_out="/tmp/ft-sort-1M.out"
rm -f ${base_line_out}
aws s3 cp s3://${AWS_BUCKET}/oneliners/inputs/1M.txt $input_file
sort $input_file -o $base_line_out

leash_out="/tmp/leash-ft-sort-1M.out"
rm -f ${leash_out}
aws s3 cp s3://${AWS_BUCKET}/ft/sort-1M.txt ${leash_out}

echo "Comparing leash output with baseline..."
if diff -q ${base_line_out} ${leash_out} > /dev/null; then
    echo "✅ Leash output matches baseline!"
else
    echo "❌ Leash output does NOT match baseline!"
fi

rm -f ${base_line_out} ${leash_out}