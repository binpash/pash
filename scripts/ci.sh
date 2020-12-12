#!/bin/bash

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

## TODO: We have to remove the sudo from the install script 
##       to be able to run it.
install() {
  scripts/install.sh
}

pash_tests() {
  cd ../compiler 
  ./test_evaluation_scripts.sh | tee  >(grep '^Summary' | cut -d ' ' -f2 > pash_tests.sum)
  PASH_RESULTS=$(cat pash_tests.sum)
  cd ../scripts
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

trap 'err_report $LINENO' ERR

echo $(date '+%F %T') "Start CI" > $RF
echo "Lots of tests"   >> $RF
pash_tests >> $RF
smoosh_tests >> $RF
echo $(date '+%F %T') "End CI"  >> $RF

RES="$(echo $PASH_RESULTS '|' $SMOOSH_RESULTS)"

FORMAT="%s  %s  %-40s  %s  %s\n"
SUM="$(printf "$FORMAT" "$(date '+%F;%T')" "$REV" "$MSG" "$RES" "$TIME")"
echo "$SUM" >> $SF


