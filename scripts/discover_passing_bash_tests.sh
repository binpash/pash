#!/bin/bash
#
# Discover which commented-out bash tests actually pass now
#

set -e

PASH_TOP="${PASH_TOP:-$(git rev-parse --show-toplevel)}"

# Use gtimeout on macOS, timeout on Linux
if command -v gtimeout &> /dev/null; then
    TIMEOUT_CMD="gtimeout"
elif command -v timeout &> /dev/null; then
    TIMEOUT_CMD="timeout"
else
    echo "Warning: timeout command not found. Install coreutils (macOS: brew install coreutils)"
    echo "Proceeding without timeouts (tests may hang)..."
    TIMEOUT_CMD=""
fi

# Helper function to run with optional timeout
run_with_timeout() {
    local seconds=$1
    shift
    if [ -n "$TIMEOUT_CMD" ]; then
        $TIMEOUT_CMD "$seconds" "$@"
    else
        "$@"
    fi
}
TEST_DIR="$PASH_TOP/evaluation/tests/interface_tests"
BASH_TESTS_DIR="$TEST_DIR/bash_tests"
RUN_SCRIPT="$TEST_DIR/run_bash_tests.sh"

PASSING_TESTS="$TEST_DIR/discovered_passing.txt"
FAILING_TESTS="$TEST_DIR/discovered_failing.txt"
RESULTS_DIR="$TEST_DIR/discovery_results"

mkdir -p "$RESULTS_DIR"
> "$PASSING_TESTS"
> "$FAILING_TESTS"

echo "Discovering passing bash tests..."
echo "================================"
echo ""

# Extract commented-out test names
COMMENTED_TESTS=$(grep -E "^    # run_test" "$RUN_SCRIPT" | sed 's/.*run_test //' | sed 's/ .*$//' | tr -d ' ')

total=0
passed=0
failed=0

for test_func in $COMMENTED_TESTS; do
    total=$((total + 1))
    
    # Extract just the test file name from the function name (test_foo.sub -> foo.sub)
    test_file=$(echo "$test_func" | sed 's/^test_//')
    
    if [ ! -f "$BASH_TESTS_DIR/$test_file" ]; then
        echo "[$total] SKIP: $test_file (file not found)"
        continue
    fi
    
    echo -n "[$total] Testing $test_file... "
    
    bash_out="$RESULTS_DIR/${test_file}.bash.out"
    pash_out="$RESULTS_DIR/${test_file}.pash.out"
    pash_err="$RESULTS_DIR/${test_file}.pash.err"
    
    # Run with bash
    cd "$BASH_TESTS_DIR"
    if run_with_timeout 30 bash "$test_file" > "$bash_out" 2>&1; then
        bash_ec=0
    else
        bash_ec=$?
    fi
    
    # Run with pash --bash
    cd "$BASH_TESTS_DIR"
    if run_with_timeout 60 "$PASH_TOP/pa.sh" --bash "$test_file" > "$pash_out" 2> "$pash_err"; then
        pash_ec=0
    else
        pash_ec=$?
    fi
    
    # Compare outputs
    if diff -q "$bash_out" "$pash_out" > /dev/null 2>&1; then
        output_match=true
    else
        output_match=false
    fi
    
    # Determine result
    if [ "$pash_ec" -eq 0 ] || [ "$pash_ec" -eq "$bash_ec" ]; then
        if [ "$output_match" = true ]; then
            echo "PASS (output matches)"
            echo "$test_func" >> "$PASSING_TESTS"
            passed=$((passed + 1))
        else
            echo "PARTIAL (compiles but output differs)"
            echo "# $test_func  # output differs" >> "$FAILING_TESTS"
            failed=$((failed + 1))
        fi
    else
        # Check if it at least compiles
        if grep -q "Parsing error" "$pash_err" 2>/dev/null; then
            echo "FAIL (parse error)"
        elif grep -q "Error" "$pash_err" 2>/dev/null; then
            echo "FAIL (runtime error)"
        else
            echo "FAIL (ec: pash=$pash_ec, bash=$bash_ec)"
        fi
        echo "$test_func" >> "$FAILING_TESTS"
        failed=$((failed + 1))
    fi
done

echo ""
echo "================================"
echo "Summary:"
echo "  Total tested: $total"
echo "  Passed: $passed"
echo "  Failed/Partial: $failed"
echo ""
echo "Results saved to:"
echo "  Passing: $PASSING_TESTS"
echo "  Failing: $FAILING_TESTS"
echo "  Outputs: $RESULTS_DIR/"
echo ""
echo "To uncomment passing tests, run:"
echo "  cat $PASSING_TESTS | while read t; do"
echo "    sed -i '' \"s/# run_test \$t/run_test \$t/\" $RUN_SCRIPT"
echo "  done"

