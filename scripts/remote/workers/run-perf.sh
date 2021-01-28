#! /bin/bash

CONTAINER_NAME=pash-perf

#docker run -itd --name pash-perf pash bin/bash

cleanup() {
  docker container stop $CONTAINER_NAME
  docker container rm $CONTAINER_NAME
}

run() {
  docker exec $CONTAINER_NAME /bin/bash -c "cd /pash/evaluation/eurosys && $@"
}

#trap 'cleanup' EXIT
trap 'echo "<<fail>>"' ERR

if [ -f lock ]; then
    echo "Busy on existing job."
else
    touch lock
    trap 'rm lock' EXIT

    run ./execute_eurosys_one_liners.sh -s
    run ./execute_unix_benchmarks.sh -s

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
