#!/usr/bin/env bash
# time: print real in seconds, to simplify parsing

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Source utils for helper functions
. "$REPO_ROOT/scripts/utils.sh"

# Always set PASH_TOP to the correct location (src/pash)
# This overrides any stale values from old PaSh installations
export PASH_TOP="$REPO_ROOT/src/pash"

bash="bash"
pash="pash -w 2"
pash_with_bash="pash -w 2 --bash"

output_dir="$REPO_ROOT/evaluation/intro/output"
rm -rf "$output_dir"
mkdir -p "$output_dir"

( cd "$REPO_ROOT/evaluation/intro/input"; ./setup.sh ) 
run_test()
{
    local test=$1
    local args="$2"
    echo -n "Running $test..."
    TIMEFORMAT="${test%%.*}:%3R" # %3U %3S"
    { time $bash "$test" > "$output_dir/$test.bash.out"; } 2> >(tee -a $output_dir/results.time_bash)
    test_bash_ec=$?
    TIMEFORMAT="%3R" # %3U %3S"
    { time $pash $args "$test" > "$output_dir/$test.pash.out"; } 2> >(tee -a $output_dir/results.time_pash)
    test_pash_ec=$?
    diff "$output_dir/$test.bash.out" "$output_dir/$test.pash.out"
    test_diff_ec=$?
    if [ $test_diff_ec -ne 0 ]; then
        echo -n "$test output mismatch"
    fi
    if [ $test_bash_ec -ne $test_pash_ec ]; then
        echo -n "$test exit code mismatch"
    fi
    if [ $test_diff_ec -ne 0 ] || [ $test_bash_ec -ne $test_pash_ec ]; then
        echo "$test are not identical" >> $output_dir/result_status
        echo -e '\t\tFAIL'
        return 1
    fi
    TIMEFORMAT="%3R"
    { time $pash_with_bash "$test" > "$output_dir/$test.pash_bash.out"; } 2> >(tee -a $output_dir/results.time_pash_bash)
    test_pash_bash_ec=$?
    diff "$output_dir/$test.bash.out" "$output_dir/$test.pash_bash.out"
    test_diff_bash_ec=$?
    if [ $test_diff_bash_ec -ne 0 ]; then
        echo -n "$test output mismatch"
    fi
    if [ $test_bash_ec -ne $test_pash_bash_ec ]; then
        echo -n "$test exit code mismatch"
    fi
    if [ $test_diff_bash_ec -ne 0 ] || [ $test_bash_ec -ne $test_pash_bash_ec ]; then
        echo "$test are not identical" >> $output_dir/result_status
        echo -e '\t\tFAIL'
        return 1
    else
        echo "$test are identical" >> $output_dir/result_status
        echo -e '\t\tOK'
        return 0
    fi
}

run_test "demo-spell.sh"
run_test "hello-world.sh"
run_test "hello-world-bash.sh" "--bash"

if type lsb_release >/dev/null 2>&1 ; then
   distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
   distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
case "$distro" in
    freebsd*)  
        # change sed to gsed
        sed () {
            gsed $@
        }
        ;;
    *)
        ;;
esac

echo "group,Bash,Pash2" > $output_dir/results.time
paste $output_dir/results.time_*  | sed 's\,\.\g' | sed 's\:\,\g' | sed 's/\t/,/' >> $output_dir/results.time

echo "Below follow the identical outputs:"
grep "are identical" "$output_dir"/result_status | awk '{print $1}'

echo "Below follow the non-identical outputs:"     
grep "are not identical" "$output_dir"/result_status | awk '{print $1}'

TOTAL_TESTS=$(cat "$output_dir"/result_status | wc -l)
PASSED_TESTS=$(grep -c "are identical" "$output_dir"/result_status)
echo "Summary: ${PASSED_TESTS}/${TOTAL_TESTS} tests passed."
