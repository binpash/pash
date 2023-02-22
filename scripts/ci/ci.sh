#!/bin/bash

## 
# This runs the majority of the core CI job, including packaging the repo
# running tests. Do not add environment installation in this script, this should
# be done manually (both for security and convenience).  No two process of this
# script can execute in parallel, nor can this process be safely interleaved
# with any process running `git` (such as pkg). The program webhook.js serves as
# a synchronization point; do not start this script if webhook.js is running on
# the same computer and accepting requests, but rather use webhook.js (as a
# daemon) to launch this script.  Otherwise you run the risk of running into
# concurrency issues.  See additional notes on webhook.js.
##

set -ex

# Placeholder for CI
REPORT_DIR=../../reports
C=5
cd ..
PASH_TOP="$PWD"
cd -

SMOOSH_RESULTS=""

trim() {
  tr -d '\n' | awk 'length > 40{$0 = substr($0, 1, 37) "..."} {print $0}'
}

build_runtime() {
  cd ../runtime 
  make
  cd $PASH_TOP/scripts
}

pash_tests() {
  cd ../evaluation/tests/ 
  ./test_evaluation_scripts.sh | tee  >(grep '^Summary' | cut -d ' ' -f2 > pash_tests.sum)
  PASH_RESULTS=$(cat pash_tests.sum)
  cd $PASH_TOP/scripts
}

smoosh_tests() {
  cd ../../smoosh
  TEST_SHELL="$PASH_TOP/pa.sh --width 2 --log_file /tmp/log_file" make -C tests veryclean
  TEST_SHELL="$PASH_TOP/pa.sh --width 2 --log_file /tmp/log_file" make -C tests  | tee >(grep 'tests passed' | cut -d ' ' -f2 > smoosh_tests.sum)
  SMOOSH_RESULTS=$(cat smoosh_tests.sum)
  cd $PASH_TOP/scripts
}

git pull

# Vars used in report summary
REV=$(git rev-parse --short HEAD)
MSG="$(git log -1 --pretty=%B | trim | head -n 1)"
RES="fail"
TIME="0s"

# Two report files
RF=$REPORT_DIR/$REV
SF=$REPORT_DIR/summary
ISF=$REPORT_DIR/summary.inv

err_report() {
	echo "Error on line $1"
  FORMAT="%s  %s  %-40s  %s  %s\n"
  SUM="$(printf "$FORMAT" "$(date '+%F;%T')" "$REV" "$MSG" "$RES" "$TIME")"
  echo "$SUM" >> $SF
}

stage() {
  echo $(date '+%F %T') $REV $1 >> $RF
}

cleanup() {
  git clean -f
}

trap 'err_report $LINENO' ERR
trap 'cleanup' EXIT

# To respect invariants of stages
mkdir -p $REPORT_DIR ../../get

echo $(date '+%F %T') $REV "Starting" > $RF
START_TIME=$(date +%s);
stage "Packaging PaSh"
./pkg.sh
stage "Building Runtime"
build_runtime >> $RF
stage "Running PaSh Tests"
pash_tests >> $RF
stage "Running Smoosh Tests"
smoosh_tests >> $RF
stage "Completing CI"
END_TIME=$(date +%s);

RES="$(echo $PASH_RESULTS '|' $SMOOSH_RESULTS)"
TIME=$(echo $((END_TIME-START_TIME)) | awk '{print int($1/60)":"int($1%60)}')

FORMAT="%s  %s  %-40s  %s  %ss\n"
SUM="$(printf "$FORMAT" "$(date '+%F;%T')" "$REV" "$MSG" "$RES" "$TIME")"
cat $SF | awk '{a[i++]=$0} END {for (j=i-1; j>=0;) print a[j--] }' > $ISF
echo "$SUM" >> $ISF
cat $ISF | awk '{a[i++]=$0} END {for (j=i-1; j>=0;) print a[j--] }' > $SF


