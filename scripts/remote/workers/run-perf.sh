#! /bin/bash

# Invariant: Docker image `pash` must exist.

CONTAINER_NAME=pash-perf

cleanup() {
  docker container stop $CONTAINER_NAME
  docker container rm $CONTAINER_NAME
}

cleanup
docker run -itd --name pash-perf pash /bin/bash

trap 'cleanup' EXIT
trap 'echo "<<fail>>"' ERR

if [ -d lock ]; then
    echo "Busy on existing job."
else
    mkdir lock
    trap 'rm -r lock' EXIT

    docker exec $CONTAINER_NAME /bin/bash -c "/pash/scripts/ci-perf.sh"

    # WARNING: If you are using the snap edition of Docker, this
    # command may break at the start of a personalized rabbit hole.
    #
    # Do not delete the dots and slashes at the end of the arguments.
    # They trigger desired merging behavior.
    # https://github.com/moby/moby/issues/31251#issuecomment-281634234
    docker cp "$CONTAINER_NAME:/tmp/results/." results/

    # The controller looks for this.
    echo '<<done>>'
fi
