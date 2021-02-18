#!/bin/bash
# Suggests envvars for use in ./make-ec2.sh

main() {
    local vpc_id="$(aws ec2 describe-vpcs --output text --query 'Vpcs[0].VpcId')";
    local key_name="$(aws ec2 describe-key-pairs --output text --query 'KeyPairs[0].KeyName')";
    local subnet="$(aws ec2 describe-subnets --output text --query 'Subnets[0].SubnetId' --filter Name=vpc-id,Values=$vpc_id)";
    local sg="$(aws ec2 describe-security-groups --output text --filter Name=vpc-id,Values=$vpc_id --query 'SecurityGroups[0].GroupId')";

    echo "export PASH_AWS_EC2_AMI='ami-0d739ceed1874f156';";
    echo "export PASH_AWS_EC2_INSTANCE_TYPE='t2.micro';";
    echo "export PASH_AWS_EC2_VPC_ID='$vpc_id';";
    echo "export PASH_AWS_EC2_KEY_NAME='$key_name';";
    echo "export PASH_AWS_EC2_SUBNET='$subnet';";
    echo "export PASH_AWS_EC2_SECURITY_GROUP='$sg';";
    echo "export PASH_AWS_EC2_DISK_SIZE_GB='10';";
}

main
