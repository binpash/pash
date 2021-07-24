./test-common.sh tr "'[a-z]' '\n'" ../tr
./test-common.sh tr "A-Z a-z" ../tr
./test-common.sh tr "-c '[A-Z]' '\n'" ../tr
./test-common.sh tr "--complement '[1-9]\n*' '[a-z][A-Z]" ../tr

./test-common.sh wc "" ../wc
./test-common.sh wc "-l" ../wc
./test-common.sh wc "-w" ../wc
./test-common.sh wc "-c" ../wc
./test-common.sh wc "-m" ../wc
./test-common.sh wc "-L" ../wc
./test-common.sh wc "-lcm" ../wc
./test-common.sh wc "-mlw" ../wc
./test-common.sh wc "-mLc" ../wc
./test-common.sh wc "-L -mc -w" ../wc
./test-common.sh wc "--bytes -c --chars -L" ../wc
./test-common.sh wc "-L --lines --words" ../wc

./test-common.sh uniq "" ../uniq
./test-common.sh uniq "-c" ../uniq 
./test-common.sh uniq "--count" ../uniq

# These tests are run during PASH_TOP/scripts/run_tests.sh
# Make sure to build the aggregators using PASH_TOP/scripts/setup-pash.sh first
#
# More tests can be added like this:
#   ./test-common.sh cmd args agg
# where 
#   cmd   - is a shell command like   uniq
#   args  - are arguements like       -c
#   agg   - is an aggregator like     ./uniq-c