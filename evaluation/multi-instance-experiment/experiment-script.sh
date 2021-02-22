#!/usr/bin/env bash

cd pash/evaluation/eurosys
nohup ../../scripts/bg-worker.sh ../../done.txt ../../output.log ./execute_eurosys_one_liners.sh -s &
