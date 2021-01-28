#! /bin/bash

# Invariant: Docker image `pash` must exist.

CONTAINER_NAME=pash-perf

cleanup() {
  docker container stop $CONTAINER_NAME
  docker container rm $CONTAINER_NAME
}

run() {
  docker exec $CONTAINER_NAME /bin/bash -c "cd /pash/evaluation/eurosys && $@"
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

    docker exec $CONTAINER_NAME /bin/bash -c "cd /pash && git pull"
    run ./execute_eurosys_one_liners.sh -m
    run ./execute_unix_benchmarks.sh -m

    ## Uncomment when ready to deal with large inputs.
    # run execute_baseline_sort.sh
    # run execute_max_temp_dish_evaluation.sh
    # run execute_web_index_dish_evaluation.sh

    # WARNING: If you are using the snap edition of Docker, this
    # command may break at the start of a personalized rabbit hole.
    #
    # Do not delete the dots and slashes at the end of the arguments.
    # They trigger desired merging behavior.
    # https://github.com/moby/moby/issues/31251#issuecomment-281634234
    docker cp "$CONTAINER_NAME:/pash/evaluation/results/." results/

    # The controller looks for this.
    echo '<<done>>'
fi
