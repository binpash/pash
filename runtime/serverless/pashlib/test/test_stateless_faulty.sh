#!/bin/bash
cd $(dirname "$0")

JOB_ID="test-ft-nfa-regex-faulty"
FOLDER_ID="ft-scripts"

cleanup() {
    trap - EXIT INT TERM
    kill -- -$$ 2>/dev/null || true
}
trap cleanup EXIT INT TERM

export LANG=C.UTF-8
export LC_ALL=C.UTF-8

aws s3 cp scripts/nfa-regex-lambda-1-faulty.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/nfa-regex-lambda-1-faulty.sh
aws s3 cp scripts/nfa-regex-lambda-2-faulty.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/nfa-regex-lambda-2-faulty.sh
aws s3 cp scripts/nfa-regex-lambda-3-faulty.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/nfa-regex-lambda-3-faulty.sh
aws s3 cp scripts/nfa-regex-lambda-4-faulty.sh s3://${AWS_BUCKET}/sls-scripts/${FOLDER_ID}/nfa-regex-lambda-4-faulty.sh
SCRIPT_ID_1="nfa-regex-lambda-1-faulty"
SCRIPT_ID_2="nfa-regex-lambda-2-faulty"
SCRIPT_ID_3="nfa-regex-lambda-3-faulty"
SCRIPT_ID_4="nfa-regex-lambda-4-faulty"

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_1}\"],
    \"is_stateless\": true,
    \"time_to_crash\": 10
    }" \
    --invocation-type Event \
    response.json

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_2}\"],
    \"is_stateless\": true
    }" \
    --invocation-type Event \
    response.json

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_3}\"],
    \"is_stateless\": true
    }" \
    --invocation-type Event \
    response.json

aws lambda invoke \
    --function-name lambda \
    --payload "{
    \"job_id\": \"${JOB_ID}\",
    \"folder_ids\": [\"${FOLDER_ID}\"],
    \"ids\": [\"${SCRIPT_ID_4}\"],
    \"is_stateless\": true
    }" \
    --invocation-type Event \
    response.json

bash scripts/nfa-regex-ec2-faulty.sh

echo "Finished execution, checking results..."

input_file="/tmp/100M.txt"
base_line_out="/tmp/ft-nfa-regex-100M.out"
baseline_out_s3_path="s3://${AWS_BUCKET}/ft/ft-nfa-regex-100M.out"
rm -f ${input_file} ${base_line_out}

# if base_line_out presents in s3, skip the generation, otherwise generate it
if aws s3 ls ${baseline_out_s3_path} > /dev/null 2>&1; then
    echo "Baseline output already exists in S3, skipping generation."
    aws s3 cp ${baseline_out_s3_path} ${base_line_out}
else
    echo "Baseline output does not exist in S3, generating it..."
    aws s3 cp s3://${AWS_BUCKET}/oneliners/inputs/100M.txt $input_file
    cat $input_file | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2' > ${base_line_out}
    aws s3 cp ${base_line_out} ${baseline_out_s3_path}
fi

leash_out="/tmp/leash-ft-nfa-regex-100M.out"
rm -f ${leash_out}
aws s3 cp s3://${AWS_BUCKET}/ft/nfa-regex-100M.txt ${leash_out}

echo "Comparing leash output with baseline..."
if diff -q ${base_line_out} ${leash_out} > /dev/null; then
    echo "✅ Leash output matches baseline!"
else
    echo "❌ Leash output does NOT match baseline!"
fi

rm -f ${base_line_out} ${leash_out}
