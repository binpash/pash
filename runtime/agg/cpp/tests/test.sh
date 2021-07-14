./test-common.sh wc "" ../wc

# These tests are run during PASH_TOP/scripts/run_tests.sh
# Make sure to build the aggregators using PASH_TOP/scripts/setup-pash.sh first
#
# More tests can be added like this:
#   ./test-common.sh cmd args agg
# where 
#   cmd   - is a shell command like   uniq
#   args  - are arguements like       -c
#   agg   - is an aggregator like     ./uniq-c