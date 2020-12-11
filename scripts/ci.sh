#!/bin/bash

set -ex

# Placeholder for CI
REPORT_DIR=../../reports
C=5

trim() {
  awk 'length > 40{$0 = substr($0, 1, 37) "..."} {print $0}' 
}

## TODO: We have to remove the sudo from the install script 
##       to be able to run it.
install() {
  scripts/install.sh
}

tests() {
  cd ../compiler 
  ./test_evaluation_scripts.sh
  cd ../
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
tests >> $RF
echo $(date '+%F %T') "End CI"  >> $RF

RES="pass"

FORMAT="%s  %s  %-40s  %s  %s\n"
SUM="$(printf "$FORMAT" "$(date '+%F;%T')" "$REV" "$MSG" "$RES" "$TIME")"
echo "$SUM" >> $SF


