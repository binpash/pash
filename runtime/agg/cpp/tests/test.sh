./test-common.sh tr "'[a-z]' '\n'" ../bin/tr
./test-common.sh tr "A-Z a-z" ../bin/tr
./test-common.sh tr "-c '[A-Z]' '\n'" ../bin/tr
./test-common.sh tr "--complement '[1-9]\n*' '[a-z][A-Z]" ../bin/tr
./test-common.sh tr "--complement -t '[1-9]\n*' '[a-z][A-Z]" ../bin/tr
./test-common.sh tr "-d '\n'" ../bin/tr
./test-common.sh tr "-tcd '[1-9][a-z][A-Z]\n'" ../bin/tr

./test-common.sh wc "" ../bin/wc
./test-common.sh wc "-l" ../bin/wc
./test-common.sh wc "-w" ../bin/wc
./test-common.sh wc "-c" ../bin/wc
./test-common.sh wc "-m" ../bin/wc
./test-common.sh wc "-L" ../bin/wc
./test-common.sh wc "-lcm" ../bin/wc
./test-common.sh wc "-mlw" ../bin/wc
./test-common.sh wc "-mLc" ../bin/wc
./test-common.sh wc "-L -mc -w" ../bin/wc
./test-common.sh wc "--bytes -c --chars -L" ../bin/wc
./test-common.sh wc "-L --lines --words" ../bin/wc

./test-common.sh uniq "" ../bin/uniq
./test-common.sh uniq "-c" ../bin/uniq 
./test-common.sh uniq "--count" ../bin/uniq

# These tests are run during PASH_TOP/scripts/run_tests.sh
# Make sure to build the aggregators using PASH_TOP/scripts/setup-pash.sh first
#
# More tests can be added like this:
#   ./test-common.sh cmd args agg
# where 
#   cmd   - is a shell command like   uniq
#   args  - are arguements like       -c
#   agg   - is an aggregator like     ./uniq-c
