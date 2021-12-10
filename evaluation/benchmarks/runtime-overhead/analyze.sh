#!/bin/bash

echo "BaSh average:"
cat total-times.txt | grep "bash:" | cut -d: -f 2 | awk "{ total += \$1; count++ } END { print total/count }"

echo "Blish -expand -server -par_pipes average:"
cat total-times.txt | grep "PaSh_no_daemon_bash_mirror:" | cut -d: -f 2 | awk "{ total += \$1; count++ } END { print total/count }"

echo "Blish -server -par_pipes average:"
cat total-times.txt | grep "PaSh_no_daemon:" | cut -d: -f 2 | awk "{ total += \$1; count++ } END { print total/count }"

echo "Blish -par_pipes average:"
cat total-times.txt | grep "PaSh_daemon:" | cut -d: -f 2 | awk "{ total += \$1; count++ } END { print total/count }"

echo "Blish average:"
cat total-times.txt | grep "PaSh_daemon_parallel_pipelines:" | cut -d: -f 2 | awk "{ total += \$1; count++ } END { print total/count }"