#!/bin/bash
# First command is almost always a generator
set -e
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/aliases/input/out}


# grep -iv ': starting\|kernel: .*: Power Button\|watching system buttons\|Stopped Cleaning Up\|Started Crash recovery kernel' /var/log/messages /var/log/syslog /var/log/* 2> /dev/null | grep -iw 'recover[a-z]*\|power[a-z]*\|shut[a-z ]*down\|rsyslogd\|ups' > /tmp/__shutdown.log && echo 'File written to /tmp__shutdown.log'
# doesn't do much :/
grep -iv ': starting\|kernel: .*: Power Button\|watching system buttons\|Stopped Cleaning Up\|Started Crash recovery kernel' /var/log/messages /var/log/syslog /var/log/* 2> /dev/null | 
    grep --regex 'recover[a-z]*\|power[a-z]*\|shut[a-z ]*down\|rsyslogd\|ups' > ${OUT}/shutdown.log 
