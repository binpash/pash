#!/usr/bin/env bash

INSTANCE_ID="${1:?No instance-id given.}"
KEY_PATH="${2:?No key path given.}"

$PASH_TOP/scripts/with-ec2.sh "${INSTANCE_ID}" $PASH_TOP/evaluation/multi-instance-experiment/ssh-script.sh "${KEY_PATH}"