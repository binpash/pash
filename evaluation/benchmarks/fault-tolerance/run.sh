#!/bin/bash
cd $(dirname "$0")

export RUST_LOG=info

stateless_in="oneliners/inputs/50M.txt"
stateless_out="oneliners/outputs/ft-stateless-"
stateless_script="oneliners/scripts/nfa-regex.sh"
stateful_in="oneliners/inputs/3G.txt"
stateful_out="oneliners/outputs/ft-stateful-"
stateful_script="oneliners/scripts/sort.sh"

fetch_logs_and_cost() {
    local logs_dir="$1"
    local start_time_ms="$2"
    local job_id="$3"

    sleep 12

    if [ -d "$logs_dir" ]; then
        echo "Removing existing logs directory: $logs_dir"
        rm -rf "$logs_dir"
    fi

    local utils_output
    utils_output=$(python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir" "$start_time_ms" "$job_id")
    echo "$utils_output"
}


run_case() {
    local case_name="$1"
    local in_path="$2"
    local out_prefix="$3"
    local script_path="$4"

    export IN="$in_path"

    echo "[$case_name] Running fault-free execution (${case_name})"
    JOB_ID="ft_${case_name}_$(date +%s%N)"
    export PASH_JOB_ID="$JOB_ID"
    export OUT="${out_prefix}baseline-"
    START_TIME_MS=$(($(date +%s%3N) - 10000))
    start_ns=$(date +%s%N)
    USE_DYNAMIC_BOUNDARIES=true NO_RESPLITTING=true PASH_S3_CHUNKS_PER_LAMBDA=16 \
        $PASH_TOP/pa.sh --serverless_exec --enable_s3_direct --no_resplitting --ec2_width $(nproc) -w4 \
        $PASH_TOP/evaluation/benchmarks/${script_path} 2>&1 >"ec2_${case_name}.log"
    end_ns=$(date +%s%N)
    faultfree_wall_clock_time=$(awk -v start="$start_ns" -v end="$end_ns" 'BEGIN{printf "%.3f", (end-start)/1000000000}')
    echo "[$case_name] Wall clock time: ${faultfree_wall_clock_time} seconds"
    sleep 15
    fetch_logs_and_cost "logs/${JOB_ID}" "$START_TIME_MS" "${JOB_ID}"

    # Goto the second log under logs/$(JOB_ID)/logs and find the lambda execution time
    lambda_time=$(grep "REPORT RequestId" logs/${JOB_ID}/logs/* | head -n 2 | tail -n 1 | awk -F "Billed Duration: " '{print $2}' | awk -F " ms" '{print $1}')
    lambda_time=$(awk -v ms="$lambda_time" 'BEGIN{printf "%.3f", ms/1000}')
    echo "[$case_name] Lambda execution time: ${lambda_time} seconds"

    for crash_progress in 20 80 ; do
        JOB_ID="ft_${case_name}_crash_${crash_progress}_$(date +%s%N)"
        export PASH_JOB_ID="$JOB_ID"
        export OUT="${out_prefix}crash-${crash_progress}-"
        START_TIME_MS=$(($(date +%s%3N) - 10000))
        start_ns=$(date +%s%N)
        export PASH_TIME_TO_CRASH=$(awk -v lt="$lambda_time" -v p="$crash_progress" 'BEGIN{printf "%d", (lt * p / 100) + 0.999999}')
        echo "[$case_name] Injecting crash at ${crash_progress}% progress (PASH_TIME_TO_CRASH=${PASH_TIME_TO_CRASH}s)"
        PASH_LAMBDA_CRASH_IDX=1 USE_DYNAMIC_BOUNDARIES=true NO_RESPLITTING=true PASH_S3_CHUNKS_PER_LAMBDA=16 \
            $PASH_TOP/pa.sh --serverless_exec --enable_s3_direct --no_resplitting --ec2_width $(nproc) -w4 \
            $PASH_TOP/evaluation/benchmarks/${script_path} 2>&1 >"ec2_${case_name}_faulty_${crash_progress}.log"
        end_ns=$(date +%s%N)
        wall_clock_time=$(awk -v start="$start_ns" -v end="$end_ns" 'BEGIN{printf "%.3f", (end-start)/1000000000}')
        echo "[$case_name] Wall clock time with crash at ${crash_progress}%: ${wall_clock_time} seconds"
        sleep 15
        fetch_logs_and_cost "logs/${JOB_ID}" "$START_TIME_MS" "${JOB_ID}"
    done
}

rm all.log
run_case "stateless" "$stateless_in" "$stateless_out" "$stateless_script" 2>&1 | tee -a all.log
run_case "stateful" "$stateful_in" "$stateful_out" "$stateful_script" 2>&1 | tee -a all.log
