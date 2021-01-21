#!/bin/bash

# Runs correctness tests and sends reports to shared storage.
#
# Invariants:
# - Docker is installed.
# - A `pash-playground` container is available as documented in pash's README.md. It doesn't matter if its running.

set -e
trap 'rm lock' EXIT
trap 'echo "<<fail>>"' ERR

if [ "$USER" == root ]; then
  runuser -l ubuntu -c $0
elif [ -f lock ]; then
  echo "Busy on existing job."
else
    touch lock
    docker start pash-playground
    docker exec pash-playground bash -c 'cd /pash/scripts && ./ci.sh'

    # The snap edition of Docker is strange when it comes to the cp command.
    #
    # First, docker cp can't write to this directory, even using sudo.
    # The workaround is to copy to /tmp.
    # https://stackoverflow.com/a/45276559
    #
    # Second, docker cp doesn't actually put it directly under /tmp.
    # It puts it under /tmp/snap.docker/tmp. If you are not using
    # the snap edition of Docker, then you may need to adjust the
    # first argument to rsync.

    sudo docker cp pash-playground:/reports /tmp
    sudo rsync -a /tmp/snap.docker/tmp/reports reports
    sudo chown -R ubuntu:ubuntu reports

    # This allows shared review of reports
    aws s3 cp --recursive reports "s3://pash-reports" \
              --acl 'public-read' \
              --content-type 'text/plain'

   echo '<<done>>'
fi
