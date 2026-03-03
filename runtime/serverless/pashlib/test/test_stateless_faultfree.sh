#!/bin/bash
cd $(dirname "$0")

JOB_ID="test-ft-nfa-regex"
FOLDER_ID="ft-scripts"

export LANG=C.UTF-8
export LC_ALL=C.UTF-8

aws s3 cp scripts/nfa-regex-lambda-1.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/nfa-regex-lambda-1.sh
aws s3 cp scripts/nfa-regex-lambda-2.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/nfa-regex-lambda-2.sh
aws s3 cp scripts/test.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/test.sh
SCRIPT_ID_1="nfa-regex-lambda-1"
SCRIPT_ID_2="nfa-regex-lambda-2"

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

bash scripts/nfa-regex-ec2.sh

echo "Finished execution, checking results..."

input_file="/tmp/1M.txt"
base_line_out="/tmp/ft-nfa-regex-1M.out"
rm -f ${base_line_out}
aws s3 cp s3://${AWS_BUCKET}/oneliners/inputs/1M.txt $input_file
cat $input_file | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > ${base_line_out}

leash_out="/tmp/leash-ft-nfa-regex-1M.out"
rm -f ${leash_out}
aws s3 cp s3://${AWS_BUCKET}/ft/nfa-regex-1M.txt ${leash_out}

echo "Comparing leash output with baseline..."
if diff -q ${base_line_out} ${leash_out} > /dev/null; then
    echo "✅ Leash output matches baseline!"
else
    echo "❌ Leash output does NOT match baseline!"
fi

rm -f ${base_line_out} ${leash_out}