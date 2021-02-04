#! /bin/bash

# Invariant: ~/pash/scripts/Dockerfile exists, if argument is not
# specified.

IMAGE_TAG=pash-perf-img
CONTAINER_NAME=pash-perf

cleanup() {
    rm -r lock
    docker container stop $CONTAINER_NAME
    docker container rm $CONTAINER_NAME
}

if [ -d lock ]; then
    echo "Busy on existing job."
else
    cleanup
    mkdir lock

    trap 'cleanup' EXIT
    trap 'echo "<<fail>>"' ERR

    docker build -t $IMAGE_TAG - < "${1:-~/pash/scripts/Dockerfile}"
    docker run -itd --name $CONTAINER_NAME $IMAGE_TAG /bin/bash

    # /!\ Remove this line when merging branch to main /!\
    docker exec $CONTAINER_NAME /bin/bash -c "cd /pash && git fetch && git checkout perf && git pull"

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
