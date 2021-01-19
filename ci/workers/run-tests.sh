#!/bin/bash

# Runs correctness tests and sends reports to shared storage.
#
# Invariants:
# - Docker is installed.
# - AWS CLI is installed, authorized to write to a pash-reports bucket on S3, and authorized to set ACLs for related objects.
# - A `pash-playground` container is available as documented in pash's README.md. It doesn't matter if its running.

trap 'rm lock' EXIT

if [ "$USER" == root ]; then
  runuser -l ubuntu -c $0
elif [ ! -f lock ]; then
    touch lock
    rm -rf reports
    docker start pash-playground
    docker exec pash-playground bash -c 'rm -rf /reports; mkdir /reports; cd /pash/scripts && ./ci.sh'
    docker cp pash-playground:/reports reports
    aws s3 cp --recursive reports "s3://pash-reports/$(date '+%F_%T' | tr : -)" \
              --acl 'public-read' \
              --content-type 'text/plain'
fi
