#!/bin/bash

# Runs correctness tests and sends reports to shared storage.
#
# Requires Docker and an available `pash-playground` container (See pash's README.md).

set -e
trap 'echo "<<fail>>"' ERR

if [ -f lock ]; then
  echo "Busy on existing job."
else
    touch lock
    trap 'rm lock' EXIT
    docker start pash-playground
    docker exec pash-playground bash -c 'cd /pash/scripts && ./ci.sh'

    # WARNING: If you are using the snap edition of Docker, this
    # command may break at the start of a personalized rabbit hole.
    #
    # Do not delete the dots and slashes at the end of the arguments.
    # They trigger desired merging behavior.
    # https://github.com/moby/moby/issues/31251#issuecomment-281634234
    docker cp pash-playground:/reports/. reports/

    # The controller looks for this.
    echo '<<done>>'
fi
