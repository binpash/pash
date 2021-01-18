#!/bin/sh

# Used by the controller to start work. In plain English, this tells
# AWS "Tell every EC2 instance tagged as a Pash worker to run their
# assigned script."
#
# The other options are just there to set safety limits and tell AWS
# more about how to do that.
#
# TODO: allow (safe) overriding of --targets.

aws ssm send-command \
    --document-name "AWS-RunShellScript" \
    --document-version "1" \
    --targets '[{"Key":"tag:Pash","Values":["worker"]}]' \
    --parameters '{"workingDirectory":["pash"],"executionTimeout":["3600"],"commands":["./worker-script.sh",""]}' \
    --timeout-seconds 600 \
    --max-concurrency "50" \
    --max-errors "0" \
    --output-s3-bucket-name "pash-reports" \
    --region us-east-1
