#!/usr/bin/env bash

DONE_FILE=${1?Done file argument not given}
KEY_PATH=${2?Key path not given}
HOSTNAME=${3?Hostname not given}

check_if_file_exists()
{
"${PASH_TOP}/evaluation/multi-instance-experiment/ssh-script.sh" "${KEY_PATH}" "${HOSTNAME}" <<EOF
if [[ -f $DONE_FILE ]]; then
    cat $DONE_FILE
    exit 0
else
    exit 1
fi
EOF
}

until check_if_file_exists
do
#   echo "Not completed"
  sleep 5
done
