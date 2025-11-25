#!/bin/bash
# Holepunch Performance Benchmarking Suite
# Tests different file sizes and output modes with statistical analysis

set -e

PASH_TOP="/home/ubuntu/pash"
ORCHESTRATOR="$PASH_TOP/minimal_orchestrate_bench.sh"

# Configuration
FILE_SIZES=("100M" "500M" "1G")
OUTPUT_MODES=("/dev/null" "/tmp/f.out")
RUNS_PER_CONFIG=9

# Results storage
declare -A ec2_times
declare -A lambda_times
declare -A total_times

# Helper: Calculate statistics (avg, min, max) from array of times
calculate_stats() {
    local times_str="$1"
    IFS=',' read -ra times <<< "$times_str"

    # Convert to milliseconds for precise arithmetic
    local times_ms=()
    for time in "${times[@]}"; do
        time_sec=${time%s}
        time_ms=$(awk "BEGIN {printf \"%.0f\", $time_sec*1000}")
        times_ms+=($time_ms)
    done

    # Calculate statistics in milliseconds
    local sum=0
    local min=${times_ms[0]}
    local max=${times_ms[0]}

    for time_ms in "${times_ms[@]}"; do
        sum=$((sum + time_ms))
        [ $time_ms -lt $min ] && min=$time_ms
        [ $time_ms -gt $max ] && max=$time_ms
    done

    local avg=$((sum / ${#times_ms[@]}))

    # Convert back to seconds for display
    local avg_sec=$(awk "BEGIN {printf \"%.3f\", $avg/1000}")
    local min_sec=$(awk "BEGIN {printf \"%.3f\", $min/1000}")
    local max_sec=$(awk "BEGIN {printf \"%.3f\", $max/1000}")

    echo "avg=${avg_sec}s  min=${min_sec}s  max=${max_sec}s"
}

# Helper: Calculate throughput in MB/s
calculate_throughput() {
    local file_size=$1
    local time_sec=$2

    # Convert file size to MB
    local size_mb
    case $file_size in
        100M) size_mb=100 ;;
        500M) size_mb=500 ;;
        1G)   size_mb=1024 ;;
    esac

    # Remove 's' suffix and calculate
    time_sec=${time_sec%s}
    local throughput=$(awk "BEGIN {printf \"%.1f\", $size_mb/$time_sec}")
    echo "${throughput}MB/s"
}

echo "================================================================="
echo "Holepunch Benchmark Suite"
echo "================================================================="
echo "File sizes: ${FILE_SIZES[*]}"
echo "Output modes: ${OUTPUT_MODES[*]}"
echo "Runs per config: $RUNS_PER_CONFIG"
echo "Total runs: $((${#FILE_SIZES[@]} * ${#OUTPUT_MODES[@]} * RUNS_PER_CONFIG))"
echo "================================================================="
echo ""

# Create results directory
mkdir -p "$PASH_TOP/benchmark_results"
RESULTS_FILE="$PASH_TOP/benchmark_results/results_$(date +%Y%m%d_%H%M%S).txt"

# Run benchmarks
config_num=0
total_configs=$((${#FILE_SIZES[@]} * ${#OUTPUT_MODES[@]}))

for file_size in "${FILE_SIZES[@]}"; do
    for output_mode in "${OUTPUT_MODES[@]}"; do
        config_num=$((config_num + 1))
        config_key="${file_size}_${output_mode//\//_}"

        echo "[$config_num/$total_configs] Running: $file_size -> $output_mode"

        ec2_times_arr=()
        total_times_arr=()

        for run in $(seq 1 $RUNS_PER_CONFIG); do
            echo -n "  Run $run/$RUNS_PER_CONFIG... "

            # Run the benchmark and capture output
            output=$("$ORCHESTRATOR" "$file_size" "$output_mode" 2>&1)
            exit_code=$?

            if [ $exit_code -eq 0 ]; then
                # Parse timing from output (Lambda timing is in CloudWatch, not available here)
                ec2_time=$(echo "$output" | grep "EC2_TIME:" | awk '{print $2}')
                total_time=$(echo "$output" | grep "TOTAL_TIME:" | awk '{print $2}')

                if [ -n "$ec2_time" ] && [ -n "$total_time" ]; then
                    ec2_times_arr+=("$ec2_time")
                    total_times_arr+=("$total_time")
                    echo "✓ (total: $total_time)"
                else
                    echo "✗ (timing data missing - got: ec2=$ec2_time total=$total_time)"
                fi
            else
                echo "✗ (exit code: $exit_code)"
            fi

            # Small delay between runs
            sleep 2
        done

        # Store results
        ec2_times[$config_key]=$(IFS=,; echo "${ec2_times_arr[*]}")
        total_times[$config_key]=$(IFS=,; echo "${total_times_arr[*]}")

        echo ""
    done
done

# Generate report
echo "================================================================="
echo "BENCHMARK RESULTS"
echo "================================================================="
echo "" | tee "$RESULTS_FILE"

for file_size in "${FILE_SIZES[@]}"; do
    for output_mode in "${OUTPUT_MODES[@]}"; do
        config_key="${file_size}_${output_mode//\//_}"

        echo "Config: $file_size -> $output_mode" | tee -a "$RESULTS_FILE"

        if [ -n "${ec2_times[$config_key]}" ]; then
            ec2_stats=$(calculate_stats "${ec2_times[$config_key]}")
            total_stats=$(calculate_stats "${total_times[$config_key]}")

            # Get min/avg/max times for throughput calculation
            IFS=',' read -ra total_arr <<< "${total_times[$config_key]}"
            min_time=$(calculate_stats "${total_times[$config_key]}" | grep -o 'min=[^ ]*' | cut -d= -f2)
            max_time=$(calculate_stats "${total_times[$config_key]}" | grep -o 'max=[^ ]*' | cut -d= -f2)
            avg_time=$(calculate_stats "${total_times[$config_key]}" | grep -o 'avg=[^ ]*' | cut -d= -f2)

            min_throughput=$(calculate_throughput "$file_size" "$max_time")
            max_throughput=$(calculate_throughput "$file_size" "$min_time")
            avg_throughput=$(calculate_throughput "$file_size" "$avg_time")

            echo "  EC2 Time:   $ec2_stats" | tee -a "$RESULTS_FILE"
            echo "  Total Time: $total_stats" | tee -a "$RESULTS_FILE"
            echo "  Throughput: avg=$avg_throughput min=$min_throughput max=$max_throughput" | tee -a "$RESULTS_FILE"
            echo "  (Lambda timing available in CloudWatch Logs)" | tee -a "$RESULTS_FILE"
        else
            echo "  No successful runs" | tee -a "$RESULTS_FILE"
        fi

        echo "" | tee -a "$RESULTS_FILE"
    done
done

echo "=================================================================" | tee -a "$RESULTS_FILE"
echo "Results saved to: $RESULTS_FILE"
echo "================================================================="
