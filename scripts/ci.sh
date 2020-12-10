#!/bin/bash

# Placeholder for CI
REPORT_DIR=../../reports
SUMMARY=$REPORT_DIR/summary.txt

run () {
REV=$(git rev-parse --short HEAD)
REPORT=$REPORT_DIR/$REV.txt
PASS="fail"
TIME="0s"
echo $(date '+%F %T') "Start CI" > $REPORT
echo "Lots of tests"   >> $REPORT
echo $(date '+%F %T') "End CI"  >> $REPORT

cat <(echo $(date '+%F %T') $REV $PASS $TIME) <(cat $SUMMARY) > $SUMMARY
}

git pull
run

