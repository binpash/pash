#!/bin/bash
# tag: distributed cat of log files
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/posh}
M0=${M0:-$PASH_TOP/evaluation/benchmarks/posh/input/cr_data/apps/distributed_logs/mount0}
M1=${M1:-$PASH_TOP/evaluation/benchmarks/posh/input/cr_data/apps/distributed_logs/mount1}
M2=${M2:-$PASH_TOP/evaluation/benchmarks/posh/input/cr_data/apps/distributed_logs/mount2}
M3=${M3:-$PASH_TOP/evaluation/benchmarks/posh/input/cr_data/apps/distributed_logs/mount3}
M4=${M4:-$PASH_TOP/evaluation/benchmarks/posh/input/cr_data/apps/distributed_logs/mount4}
cat ${M0}/*.csv ${M1}/*.csv ${M2}/*.csv ${M3}/*.csv ${M4}/*.csv | grep '128.151.150' 
