#! /bin/bash

# Runs a command on the host within the lifespan of an EC2 instance.
# Uses AWS CLI v2

get-instance-ip() {
    aws ec2 describe-instances \
        --instance-ids "$1" \
        --query 'Reservations[0].Instances[0].PublicIpAddress' \
        --output text
}

# This is necessary because a public IP might not be
# available immediately after the instance starts.
wait-for-instance-ip() {
    local variant="$(get-instance-ip "$1")";

    if [[ "$variant" =~ [0-9] ]]; then
        echo "$variant";
    else
        sleep 3;
        wait-for-instance-ip "$1"
    fi
}

call-with-active-ec2() {
    set -e
    local instance_id="$1"
    aws ec2 start-instances --instance-ids "$instance_id"
    local ip=$(wait-for-instance-ip "$instance_id");
    echo "$ip"
    trap "aws ec2 stop-instances --instance-ids $instance_id" EXIT
    ${@:2} "$ip"
}

# First expression detects if script is being sourced.
# https://stackoverflow.com/a/28776166
(return 0 2>/dev/null) || call-with-active-ec2 "$@"
