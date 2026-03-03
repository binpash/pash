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

echo "Finished execution, checking results..."

input_file="/tmp/1G.txt"
base_line_out="/tmp/ft-sort-1G.out"
baseline_out_s3_path="s3://${AWS_BUCKET}/ft/ft-sort-1G.out"
rm -f ${input_file} ${base_line_out}

# if base_line_out presents in s3, skip the generation, otherwise generate it
if aws s3 ls ${baseline_out_s3_path} > /dev/null 2>&1; then
    echo "Baseline output already exists in S3, skipping generation."
    aws s3 cp ${baseline_out_s3_path} ${base_line_out}
else
    echo "Baseline output does not exist in S3, generating it..."
    aws s3 cp s3://${AWS_BUCKET}/oneliners/inputs/1G.txt $input_file
    sort $input_file -o $base_line_out
    aws s3 cp ${base_line_out} ${baseline_out_s3_path}
fi

leash_out="/tmp/leash-ft-sort-1G.out"
rm -f ${leash_out}
aws s3 cp s3://${AWS_BUCKET}/ft/1G.txt ${leash_out}

echo "Comparing leash output with baseline..."
if diff -q ${base_line_out} ${leash_out} > /dev/null; then
    echo "✅ Leash output matches baseline!"
else
    echo "❌ Leash output does NOT match baseline!"
fi

rm -f ${base_line_out} ${leash_out}
