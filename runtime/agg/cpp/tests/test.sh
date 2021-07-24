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
./test-common.sh uniq " -c" ../uniq       # Works with both short
./test-common.sh uniq " --count" ../uniq  # and long option

# These tests are run during PASH_TOP/scripts/run_tests.sh
# Make sure to build the aggregators using PASH_TOP/scripts/setup-pash.sh first
#
# More tests can be added like this:
#   ./test-common.sh cmd args agg
# where 
#   cmd   - is a shell command like   uniq
#   args  - are arguements like       -c
#   agg   - is an aggregator like     ./uniq-c