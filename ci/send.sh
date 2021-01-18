#!/bin/sh

# This distributes testing work to EC2 instances tagged as test
# workers.  Run using the webhooks controller.

aws ssm send-command \
    --document-name "AWS-RunShellScript" \
    --document-version "1" \
    --targets '[{"Key":"tag:Pash","Values":["pash-test-worker"]}]' \
    --parameters '{"workingDirectory":["pash"],"executionTimeout":["3600"],"commands":["./scripts/ci.sh",""]}' \
    --timeout-seconds 600 \
    --max-concurrency "50" \
    --max-errors "0" \
    --output-s3-bucket-name "pash-reports" \
    --region us-east-1
