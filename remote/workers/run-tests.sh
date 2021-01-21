#!/bin/bash

# Runs correctness tests and sends reports to shared storage.
#
# Invariants:
# - Docker is installed.
# - AWS CLI is installed, authorized to write to a pash-reports bucket on S3, and authorized to set ACLs for related objects.
# - A `pash-playground` container is available as documented in pash's README.md. It doesn't matter if its running.

set -e
trap 'rm lock' EXIT

if [ "$USER" == root ]; then
  runuser -l ubuntu -c $0
elif [ ! -f lock ]; then
    touch lock
    docker start pash-playground
    docker exec pash-playground bash -c 'cd /pash/scripts && ./ci.sh'

    # The snap edition of Docker is strange when it comes to the cp command.
    #
    # First, you can't run docker cp into this directory, even using sudo.
    # The workaround is to copy to /tmp.
    # https://stackoverflow.com/a/45276559
    #
    # Second, docker cp doesn't actually put it directly under /tmp.
    # It puts it under /tmp/snap.docker/tmp. I have no explanation as to why.

    sudo docker cp pash-playground:/reports /tmp
    sudo cp -r /tmp/snap.docker/tmp/reports reports
    sudo chown -R ubuntu:ubuntu reports

    aws s3 cp --recursive reports "s3://pash-reports" \
              --acl 'public-read' \
              --content-type 'text/plain'
fi
