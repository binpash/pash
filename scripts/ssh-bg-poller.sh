#!/usr/bin/env bash

DONE_FILE=${1?Done file argument not given}
LOG_FILE=${2?Done file argument not given}
KEY_PATH=${3?Key path not given}
HOSTNAME=${4?Hostname not given}

check_if_file_exists()
{
"${PASH_TOP}/evaluation/multi-instance-experiment/ssh-script.sh" "${KEY_PATH}" "${HOSTNAME}" <<EOF
if [[ -f $DONE_FILE ]]; then
    cat $DONE_FILE
    exit 0
else
    cat $LOG_FILE
    exit 1
fi
EOF
}

LOCAL_LOG="$(mktemp -u)"
PREVIOUS_LOCAL_LOG="$(mktemp -u)"
touch $PREVIOUS_LOCAL_LOG

until check_if_file_exists > $LOCAL_LOG
do
  comm -23 --nocheck-order $LOCAL_LOG $PREVIOUS_LOCAL_LOG
  cat $LOCAL_LOG > $PREVIOUS_LOCAL_LOG
  sleep 5
done
