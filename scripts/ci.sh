#!/bin/bash

# Placeholder for CI
REPORT_DIR=../../reports
C=5

trim() {
  awk 'length > 40{$0 = substr($0, 1, 37) "..."} {printf "%-20s", $0}' 
}

run () {
  git pull

  #all this runs after pull
  REV=$(git rev-parse --short HEAD)
  MSG="$(git log -1 --pretty=%B | trim)"
  sleep 30

  RF=$REPORT_DIR/$REV
  SF=$REPORT_DIR/summary
  PASS="fail"
  TIME="0s"
  echo $(date '+%F %T') "Start CI" > $RF
  echo "Lots of tests"   >> $RF
  echo $(date '+%F %T') "End CI"  >> $RF

  FORMAT="%s  %s  %s  %s  %s\n"
  SUM="$(printf "$FORMAT" "$(date '+%F;%T')" "$REV" "$MSG" "$PASS" "$TIME")"
  cat <(echo "$SUM") <(cat $SF) > $SF
}

run

