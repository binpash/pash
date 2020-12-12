#!/bin/bash
## Talk to [Nikos](nikos@vasilak.is) if you need to change this script.
## Do not start this script if webhook.js is running on the same server, 
## Otherwise you run the risk of running into concurrency issues.

set -ex

# Placeholder for CI
REPORT_DIR=../../reports
C=5
cd ..
PASH_TOP="$PWD"
cd -

SMOOSH_RESULTS=""

trim() {
  awk 'length > 40{$0 = substr($0, 1, 37) "..."} {print $0}' 
}

build_runtime() {
  cd ../runtime 
  make
  cd $PASH_TOP/scripts
}

pash_tests() {
  cd ../compiler 
  ./test_evaluation_scripts.sh | tee  >(grep '^Summary' | cut -d ' ' -f2 > pash_tests.sum)
  PASH_RESULTS=$(cat pash_tests.sum)
  cd $PASH_TOP/scripts
}

smoosh_tests() {
  cd ../../smoosh
  TEST_SHELL="python3.8 $PASH_TOP/compiler/pash.py --split_fan_out 2 --log_file /tmp/log_file" make -C tests veryclean
  TEST_SHELL="python3.8 $PASH_TOP/compiler/pash.py --split_fan_out 2 --log_file /tmp/log_file" make -C tests  | tee >(grep 'tests passed' | cut -d ' ' -f2 > smoosh_tests.sum)
  SMOOSH_RESULTS=$(cat smoosh_tests.sum)
  cd $PASH_TOP/scripts
}

git pull

# Vars used in report summary
REV=$(git rev-parse --short HEAD)
MSG="$(git log -1 --pretty=%B | trim)"
RES="fail"
TIME="0s"

# Two report files
RF=$REPORT_DIR/$REV
SF=$REPORT_DIR/summary

err_report() {
	echo "Error on line $1"
  FORMAT="%s  %s  %-40s  %s  %s\n"
  SUM="$(printf "$FORMAT" "$(date '+%F;%T')" "$REV" "$MSG" "$RES" "$TIME")"
  echo "$SUM" >> $SF
}

stage() {
  echo $(date '+%F %T') $REV $1 >> $RF
}

trap 'err_report $LINENO' ERR

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
echo "$SUM" >> $SF


