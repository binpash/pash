#! /bin/bash

# Pair with ./suggest-ec2.sh

main() {
    set -x
    aws ec2 run-instances \
        --output text \
        --query "Instances[0].InstanceId" \
        --image-id "$PASH_AWS_EC2_AMI" \
        --instance-type "$PASH_AWS_EC2_INSTANCE_TYPE" \
        --key-name "$PASH_AWS_EC2_KEY_NAME" \
        --security-group-ids "$PASH_AWS_EC2_SECURITY_GROUP" \
        --monitoring "Enabled=false" \
        --subnet-id "$PASH_AWS_EC2_SUBNET" \
        --query 'Instances[0].InstanceId' \
        --output text
}

main
