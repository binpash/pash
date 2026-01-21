#!/bin/bash

cd "$(dirname "$0")" || exit 1


if [[ "$*" == *"--small"* ]]
then
    SCRIPT_INPUT_WIDTH=(
        # "1.sh:1_1M.txt"
        # # "2.sh:1_1M.txt"
        # "3.sh:1_1M.txt"
        # "4.sh:1_1M.txt"
        # "5.sh:2_1M.txt"
        "6.sh:3_1M.txt"
        # "7.sh:4_1M.txt"
        # "8.sh:4_1M.txt"
        # "9.sh:4_1M.txt"
        # "10.sh:4_1M.txt"
        # # "11.sh:4_1M.txt"
        # "13.sh:5_1M.txt"
        # # "14.sh:6_1M.txt"
        # "15.sh:7_1M.txt"
        # "17.sh:7_1M.txt"
        # "18.sh:8_1M.txt"
        # "19.sh:8_1M.txt"
        # "20.sh:8_1M.txt"
        # # "21.sh:8_1M.txt"
        # "23.sh:9.1_1M.txt"
        # "24.sh:9.2_1M.txt"
        # "25.sh:9.3_1M.txt"
        # "26.sh:9.4_1M.txt"
        # "28.sh:9.6_1M.txt"
        # "29.sh:9.7_1M.txt"
        # "30.sh:9.8_1M.txt"
        # "31.sh:9.9_1M.txt"
        # "32.sh:10_1M.txt"
        # "33.sh:10_1M.txt"
        # "35.sh:11_1M.txt"
    )
    INPUT_TYPE=".small"
elif [[ "$*" == *"--large"* ]]; then
    SCRIPT_INPUT_WIDTH=(
        # "1.sh:1_20G.txt"
        # "3.sh:1_20G.txt" # Error
        # "2.sh:1_20G.txt" # 4096 Mem usage
        # "4.sh:1_20G.txt" # 4096 Mem usage
        # "5.sh:2_20G.txt"
        "6.sh:3_20G.txt"
        # "7.sh:4_20G.txt"
        # "8.sh:4_20G.txt"
        # "9.sh:4_20G.txt"
        # "10.sh:4_20G.txt" # sort

        # "11.sh:4_20G.txt" # sort
        # 12.sh:4_20G.txt:16
        # "13.sh:5_20G.txt"
        # "14.sh:6_20G.txt" # sort
        # "15.sh:7_20G.txt" 
        # "17.sh:7_20G.txt" # sort
        # "18.sh:8_20G.txt"
        # "19.sh:8_20G.txt"
        # # "20.sh:8_20G.txt"
        # "21.sh:8_20G.txt"
        # "23.sh:9.1_20G.txt:32"
        # "24.sh:9.2_20G.txt"
        # "25.sh:9.3_20G.txt"
        # "26.sh:9.4_20G.txt"
        #"28.sh:9.6_20G.txt:16"
        # "29.sh:9.7_20G.txt"
        # "30.sh:9.8_20G.txt"
        # "31.sh:9.9_20G.txt"
        # "32.sh:10_20G.txt"
        # "33.sh:10_20G.txt"
        # "35.sh:11_20G.txt"

        # "34.sh:10_20G.txt"
    )
    INPUT_TYPE=".large"
    else 
       SCRIPT_INPUT_WIDTH=(
        # "1.sh:1_1G.txt"
        # "3.sh:1_1G.txt" # Error
        # "2.sh:1_1G.txt" # 4096 Mem usage
        # "4.sh:1_1G.txt" # 4096 Mem usage
        # "5.sh:2_1G.txt"
        "6.sh:3_1G.txt"
        # "7.sh:4_1G.txt"
        # "8.sh:4_1G.txt"
        # "9.sh:4_1G.txt"
        # "10.sh:4_1G.txt" # sort

        # "11.sh:4_1G.txt" # sort
        # 12.sh:4_1G.txt:16
        # "13.sh:5_1G.txt"
        # "14.sh:6_1G.txt" # sort
        # "15.sh:7_1G.txt" 
        # "17.sh:7_1G.txt" # sort
        # "18.sh:8_1G.txt"
        # "19.sh:8_1G.txt"
        # # "20.sh:8_1G.txt"
        # "21.sh:8_1G.txt"
        # "23.sh:9.1_1G.txt:32"
        # "24.sh:9.2_1G.txt"
        # "25.sh:9.3_1G.txt"
        # "26.sh:9.4_1G.txt"
        #"28.sh:9.6_1G.txt:16"
        # "29.sh:9.7_1G.txt"
        # "30.sh:9.8_1G.txt"
        # "31.sh:9.9_1G.txt"
        # "32.sh:10_1G.txt"
        # "33.sh:10_1G.txt"
        # "35.sh:11_1G.txt"

        # "34.sh:10_1G.txt"
    )
    INPUT_TYPE=""
       
fi

# Auto-detect the benchmark directory name
BENCHMARK_DIR=$(basename "$PWD")

# Check AWS_BUCKET is set
if [ -z "$AWS_BUCKET" ]; then
    echo "Error: AWS_BUCKET environment variable is not set"
    exit 1
fi

# SCRIPT_INPUT_WIDTH=(
#     #"1.sh:in_tiny.csv:16"
#     #"2.sh:in_tiny.csv:16"
#     #"3.sh:in_tiny.csv:16" # fails with addr in use err
#     # "4.sh:in_tiny.csv:16"
    
# )

# Helper function to download S3 output
download_s3_output() {
    local s3_key=$1
    local local_file=$2

    echo "Downloading s3://$AWS_BUCKET/$s3_key to $local_file"

    if aws s3 cp "s3://$AWS_BUCKET/$s3_key" "$local_file" --no-progress >/dev/null; then
        echo "✓ Successfully downloaded $s3_key"
        return 0
    else
        echo "✗ Failed to download $s3_key"
        return 1
    fi
}

# Helper function to compare two output files
compare_outputs() {
    local file1=$1
    local file2=$2
    local label1=$3
    local label2=$4

    if [ ! -f "$file1" ]; then
        echo "✗ COMPARISON FAILED: $file1 does not exist"
        return 1
    fi

    if [ ! -f "$file2" ]; then
        echo "✗ COMPARISON FAILED: $file2 does not exist"
        return 1
    fi

    echo "Comparing outputs:"
    echo "  - $label1: $file1 ($(wc -l < "$file1") lines)"
    echo "  - $label2: $file2 ($(wc -l < "$file2") lines)"
    echo ""

    if diff -q "$file1" "$file2" > /dev/null 2>&1; then
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "✓✓✓ OUTPUTS MATCH: Results are identical ✓✓✓"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        return 0
    else
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "✗✗✗ OUTPUTS DIFFER: Results are NOT identical ✗✗✗"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""
        echo "Full diff output:"
        echo "────────────────────────────────────────────────────────────────────────"
        diff -u "$file1" "$file2" || true
        echo "────────────────────────────────────────────────────────────────────────"
        echo ""
        local diff_lines=$(diff "$file1" "$file2" | wc -l)
        echo "Total diff lines: $diff_lines"
        return 1
    fi
}

# Array to track comparison results
comparison_results=()

# Run benchmarks with both optimization modes
for SCRIPT_INPUT in "${SCRIPT_INPUT_WIDTH[@]}"; do
    echo "========================================================================"
    echo "Running benchmark for $SCRIPT_INPUT"
    SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
    INPUT=$(echo "$SCRIPT_INPUT" | cut -d: -f2)
    WIDTH=$(echo "$SCRIPT_INPUT" | cut -d: -f3)
    WIDTH=${WIDTH:-16}

    # Mode 1: WITHOUT optimization (baseline)
    echo "------------------------------------------------------------------------"
    echo "[MODE 1/2] WITHOUT S3 direct streaming optimization"
    echo "Running $SCRIPT with input $INPUT and width $WIDTH"
    echo "------------------------------------------------------------------------"

    # Capture start timestamp (milliseconds since epoch)
    # Subtract 10 seconds (10000ms) to account for clock skew between local machine and AWS
    START_TIME_MS=$(($(date +%s%3N) - 10000))
    echo "Start timestamp: $START_TIME_MS (with 10s safety buffer for clock skew)"

    # time USE_SMART_BOUNDARIES=false IN="$BENCHMARK_DIR/inputs/$INPUT" OUT="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:s3opt" \
    #     $PASH_TOP/pa.sh --serverless_exec --enable_s3_direct -w"$WIDTH" scripts/"$SCRIPT"

    # exit

    time  IN="$BENCHMARK_DIR/inputs/$INPUT" OUT="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:noopt" \
        $PASH_TOP/pa.sh --serverless_exec  -w"$WIDTH" scripts/"$SCRIPT" #\
        #--graphviz svg --graphviz_dir "$PWD/pash_graphviz_noopt_$(date +%y-%m-%d-%H:%M:%S)"

    
    sleep 20

    logs_dir="logs/$SCRIPT:$INPUT:$WIDTH:noopt"
    if [ -d "$logs_dir" ]; then
        echo "Removing existing logs directory: $logs_dir"
        rm -rf "$logs_dir"
    fi
    python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir" "$START_TIME_MS"

    # Download noopt output from S3
    echo ""
    echo "Downloading noopt output from S3..."
    noopt_s3_key="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:nooptstdout.txt"
    noopt_local_file="/tmp/compare_noopt_${SCRIPT//\//_}_${INPUT}_${WIDTH}.txt"
    download_s3_output "$noopt_s3_key" "$noopt_local_file"

    

    echo ""
    echo "------------------------------------------------------------------------"
    echo "[MODE 2.1/2] WITH S3 direct streaming optimization (--enable_s3_direct) and smart boundaries"
    echo "Running $SCRIPT with input $INPUT and width $WIDTH"
    echo "------------------------------------------------------------------------"

    # Capture start timestamp (milliseconds since epoch)
    # Subtract 10 seconds (10000ms) to account for clock skew between local machine and AWS
    START_TIME_MS=$(($(date +%s%3N) - 10000))
    echo "Start timestamp: $START_TIME_MS (with 10s safety buffer for clock skew)"

    time USE_SMART_BOUNDARIES=true IN="$BENCHMARK_DIR/inputs/$INPUT" OUT="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:s3opt" \
        $PASH_TOP/pa.sh --serverless_exec --enable_s3_direct -w"$WIDTH" scripts/"$SCRIPT" #\
        #--graphviz svg --graphviz_dir "$PWD/pash_graphviz_s3opt_$(date +%y-%m-%d-%H:%M:%S)"

    sleep 20

    logs_dir="logs/$SCRIPT:$INPUT:$WIDTH:s3opt"
    if [ -d "$logs_dir" ]; then
        echo "Removing existing logs directory: $logs_dir"
        rm -rf "$logs_dir"
    fi
    python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir" "$START_TIME_MS"

    # Download s3opt output from S3
    echo ""
    echo "Downloading s3opt output from S3..."
    s3opt_s3_key="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:s3optstdout.txt"
    s3opt_local_file="/tmp/compare_s3opt_${SCRIPT//\//_}_${INPUT}_${WIDTH}.txt"
    download_s3_output "$s3opt_s3_key" "$s3opt_local_file"


    echo ""
    echo "------------------------------------------------------------------------"
    echo "[MODE 2.2/2] WITH S3 direct streaming optimization (--enable_s3_direct) and approximate non smart boundaries"
    echo "Running $SCRIPT with input $INPUT and width $WIDTH"
    echo "------------------------------------------------------------------------"

    # Capture start timestamp (milliseconds since epoch)
    # Subtract 10 seconds (10000ms) to account for clock skew between local machine and AWS
    START_TIME_MS=$(($(date +%s%3N) - 10000))
    echo "Start timestamp: $START_TIME_MS (with 10s safety buffer for clock skew)"

    time USE_SMART_BOUNDARIES=false IN="$BENCHMARK_DIR/inputs/$INPUT" OUT="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:s3optnonsmart" \
        $PASH_TOP/pa.sh --serverless_exec --enable_s3_direct -w"$WIDTH" scripts/"$SCRIPT" #\
        #--graphviz svg --graphviz_dir "$PWD/pash_graphviz_s3opt_$(date +%y-%m-%d-%H:%M:%S)"

    sleep 20

    logs_dir="logs/$SCRIPT:$INPUT:$WIDTH:s3optnonsmart"
    if [ -d "$logs_dir" ]; then
        echo "Removing existing logs directory: $logs_dir"
        rm -rf "$logs_dir"
    fi
    python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir" "$START_TIME_MS"

    # Download s3optnonsmart output from S3
    echo ""
    echo "Downloading s3optnonsmart output from S3..."
    s3optnonsmart_s3_key="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:s3optnonsmartstdout.txt"
    s3optnonsmart_local_file="/tmp/compare_s3optnonsmart_${SCRIPT//\//_}_${INPUT}_${WIDTH}.txt"
    download_s3_output "$s3optnonsmart_s3_key" "$s3optnonsmart_local_file"

    # Compare outputs
    echo ""
    echo "╔════════════════════════════════════════════════════════════════════════╗"
    echo "║  LIVE COMPARISON FOR: $SCRIPT_INPUT (noopt vs s3opt)"
    echo "╚════════════════════════════════════════════════════════════════════════╝"
    echo ""
    if compare_outputs "$noopt_local_file" "$s3opt_local_file" "noopt" "s3opt"; then
        comparison_results+=("$SCRIPT_INPUT (smart): ✓ MATCH")
        echo ""
        echo "╔════════════════════════════════════════════════════════════════════════╗"
        echo "║  ✓ RESULT: OUTPUTS MATCH for $SCRIPT_INPUT (noopt vs s3opt)"
        echo "╚════════════════════════════════════════════════════════════════════════╝"
    else
        comparison_results+=("$SCRIPT_INPUT (smart): ✗ DIFFER")
        echo ""
        echo "╔════════════════════════════════════════════════════════════════════════╗"
        echo "║  ✗ RESULT: OUTPUTS DIFFER for $SCRIPT_INPUT (noopt vs s3opt)"
        echo "╚════════════════════════════════════════════════════════════════════════╝"
    fi

    # Compare noopt with s3optnonsmart
    echo ""
    echo "╔════════════════════════════════════════════════════════════════════════╗"
    echo "║  LIVE COMPARISON FOR: $SCRIPT_INPUT (noopt vs s3optnonsmart)"
    echo "╚════════════════════════════════════════════════════════════════════════╝"
    echo ""
    if compare_outputs "$noopt_local_file" "$s3optnonsmart_local_file" "noopt" "s3optnonsmart"; then
        comparison_results+=("$SCRIPT_INPUT (nonsmart): ✓ MATCH")
        echo ""
        echo "╔════════════════════════════════════════════════════════════════════════╗"
        echo "║  ✓ RESULT: OUTPUTS MATCH for $SCRIPT_INPUT (noopt vs s3optnonsmart)"
        echo "╚════════════════════════════════════════════════════════════════════════╝"
    else
        comparison_results+=("$SCRIPT_INPUT (nonsmart): ✗ DIFFER")
        echo ""
        echo "╔════════════════════════════════════════════════════════════════════════╗"
        echo "║  ✗ RESULT: OUTPUTS DIFFER for $SCRIPT_INPUT (noopt vs s3optnonsmart)"
        echo "╚════════════════════════════════════════════════════════════════════════╝"
    fi

    echo ""
    echo "========================================================================"
    echo "Completed benchmark for $SCRIPT_INPUT (all modes)"
    echo "========================================================================"
    echo ""
done

echo ""
echo "========================================================================"
echo "ALL BENCHMARKS COMPLETED"
echo "========================================================================"
echo "Results:"
echo "  - Baseline (no opt): $BENCHMARK_DIR/outputs/*:noopt and logs/*:noopt"
echo "  - Optimized:         $BENCHMARK_DIR/outputs/*:s3opt and logs/*:s3opt"
echo "  - Graphviz graphs:   pash_graphviz_noopt_* and pash_graphviz_s3opt_*"
echo "========================================================================"
echo ""
echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                      FINAL COMPARISON SUMMARY                          ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

# Count matches and differences
match_count=0
differ_count=0
for result in "${comparison_results[@]}"; do
    if [[ "$result" == *"✓ MATCH"* ]]; then
        ((match_count++))
        echo "  ✓ $result"
    else
        ((differ_count++))
        echo "  ✗ $result"
    fi
done

echo ""
echo "────────────────────────────────────────────────────────────────────────"
echo "Total: ${#comparison_results[@]} benchmarks | ✓ Matches: $match_count | ✗ Differs: $differ_count"
echo "────────────────────────────────────────────────────────────────────────"

if [ $differ_count -eq 0 ]; then
    echo ""
    echo "╔════════════════════════════════════════════════════════════════════════╗"
    echo "║  🎉 SUCCESS: All outputs match! Optimization is correct. 🎉            ║"
    echo "╚════════════════════════════════════════════════════════════════════════╝"
else
    echo ""
    echo "╔════════════════════════════════════════════════════════════════════════╗"
    echo "║  ⚠️  WARNING: Some outputs differ! Review differences above. ⚠️         ║"
    echo "╚════════════════════════════════════════════════════════════════════════╝"
fi

# Clean up temporary comparison files
echo ""
echo "Cleaning up temporary files..."
rm -f /tmp/compare_noopt_*.txt /tmp/compare_s3opt_*.txt /tmp/compare_s3optnonsmart_*.txt
echo "Done!"
