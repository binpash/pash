#!/usr/bin/env bash

cd pash/evaluation/eurosys
## WARNING: It is necessary to close file descriptors so that ssh can exit
nohup ../../scripts/bg-worker.sh ../../done.txt ../../output.log ./execute_eurosys_one_liners.sh -s </dev/null >/dev/null 2>/dev/null &
