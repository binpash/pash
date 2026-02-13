#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_NAME=""
RUNNER_ARGS=()

while (($#)); do
    case "$1" in
        --benchmark)
            if (($# < 2)); then
                echo "Error: --benchmark requires a value" >&2
                exit 2
            fi
            BENCHMARK_NAME="$2"
            shift 2
            ;;
        --benchmark=*)
            BENCHMARK_NAME="${1#*=}"
            shift
            ;;
        --help|-h)
            cat <<'EOF'
Usage:
  run-leash-benchmark-matrix.sh --benchmark <benchmark-folder> [runner flags...]
  run-leash-benchmark-matrix.sh <benchmark-folder> [runner flags...]

Runner flags include:
  --noopt, --smart-prealigned, --approx-tail,
  --approx-dynamic, --approx-adaptive-gap, --approx-adaptive-simple,
  --approx-adaptive-single-shot, --small/--medium/--large, --skip-logs, --debug, --repeats N
  Short approx aliases: --approx-dyn, --approx-gap, --approx-simple, --approx-single, --approx-ss

Design:
  - This script is the only runner implementation.
  - Each benchmark folder only defines SCRIPT_INPUT_WIDTH selection in:
      <benchmark>/leash-matrix-inputs.sh
EOF
            exit 0
            ;;
        *)
            if [[ -z "${BENCHMARK_NAME}" ]]; then
                if [[ "$1" == --* ]]; then
                    echo "Error: mode not supported: $1" >&2
                    exit 2
                fi
                BENCHMARK_NAME="$1"
            else
                RUNNER_ARGS+=("$1")
            fi
            shift
            ;;
    esac
done

if [[ -z "${BENCHMARK_NAME}" ]]; then
    echo "Error: benchmark folder is required (e.g. --benchmark unix50)" >&2
    exit 2
fi

BENCHMARK_PATH="${SCRIPT_DIR}/${BENCHMARK_NAME}"
BENCHMARK_INPUT_CONFIG="${BENCHMARK_PATH}/leash-matrix-inputs.sh"

if [[ ! -d "${BENCHMARK_PATH}" ]]; then
    echo "Error: benchmark directory not found: ${BENCHMARK_PATH}" >&2
    exit 2
fi

if [[ ! -f "${BENCHMARK_INPUT_CONFIG}" ]]; then
    echo "Error: benchmark input config not found: ${BENCHMARK_INPUT_CONFIG}" >&2
    exit 2
fi

normalize_runner_mode_aliases() {
    local i
    local arg
    for i in "${!RUNNER_ARGS[@]}"; do
        arg="${RUNNER_ARGS[$i]}"
        case "$arg" in
            --approx-dyn|dynamic)
                RUNNER_ARGS[$i]="--approx-dynamic"
                ;;
            --approx-gap|gap)
                RUNNER_ARGS[$i]="--approx-adaptive-gap"
                ;;
            --approx-simple|simple)
                RUNNER_ARGS[$i]="--approx-adaptive-simple"
                ;;
            --approx-single|--approx-ss|single-shot)
                RUNNER_ARGS[$i]="--approx-adaptive-single-shot"
                ;;
        esac
    done
}

cd "${BENCHMARK_PATH}" || exit 1
normalize_runner_mode_aliases
set -- "${RUNNER_ARGS[@]}"

MODE_FLAGS=(
    --noopt
    --smart-prealigned
    --approx-tail
    --approx-dynamic
    --approx-adaptive-gap
    --approx-adaptive-simple
    --approx-adaptive-single-shot
)

ALL_ALLOWED_FLAGS=(
    "${MODE_FLAGS[@]}"
    --small
    --medium
    --large
    --skip-logs
    --debug
    --repeats
)

is_allowed_runner_flag() {
    local candidate="$1"
    for allowed in "${ALL_ALLOWED_FLAGS[@]}"; do
        if [[ "$candidate" == "$allowed" ]]; then
            return 0
        fi
    done
    return 1
}

validate_runner_flags() {
    local i=0
    local arg
    local repeats_value

    while [ "$i" -lt "${#RUNNER_ARGS[@]}" ]; do
        arg="${RUNNER_ARGS[$i]}"

        if [[ "$arg" == "--repeats" ]]; then
            if [ $((i + 1)) -ge "${#RUNNER_ARGS[@]}" ]; then
                echo "Error: --repeats requires a numeric value" >&2
                exit 2
            fi
            repeats_value="${RUNNER_ARGS[$((i + 1))]}"
            if ! [[ "$repeats_value" =~ ^[0-9]+$ ]]; then
                echo "Error: --repeats requires a numeric value, got: $repeats_value" >&2
                exit 2
            fi
            i=$((i + 2))
            continue
        fi

        if [[ "$arg" == --* ]] && ! is_allowed_runner_flag "$arg"; then
            echo "Error: mode not supported: $arg" >&2
            echo "Supported mode flags: ${MODE_FLAGS[*]}" >&2
            exit 2
        fi

        i=$((i + 1))
    done
}

validate_runner_flags

# Usage:
#   --noopt                  : Baseline (EC2 split + Lambda compute)
#   --smart-prealigned       : Smart prealigned chunks (EC2 scans boundaries)
#   --approx-tail            : Approx chunks + Lambda tail coordination (legacy)
#   --approx-dynamic         : Approx chunks + dynamic correction window
#   --approx-adaptive-gap    : Approx chunks + adaptive gap-window (EC2-side)
#   --approx-adaptive-simple : Approx chunks + fixed window from simple sampling
#   --approx-adaptive-single-shot : Approx chunks + adaptive single midpoint sample window
#   --dynamic/--gap/--simple/--single-shot : Short aliases for approx modes
#   --skip-logs              : Skip CloudWatch log fetching and reduce sleep time
#   --debug                  : Enable PASH debug output
#   --small/medium/large or nothing : Input size selection
#   --repeats N              : Run each non-baseline mode N times (default: 1)

# Parse mode flags
RUN_NOOPT=false
RUN_SMART_PREALIGNED=false
RUN_APPROX_TAIL=false
RUN_APPROX_DYNAMIC=false
RUN_APPROX_ADAPTIVE_GAP=false
RUN_APPROX_ADAPTIVE_SIMPLE=false
RUN_APPROX_ADAPTIVE_SINGLE_SHOT=false
PASH_DEBUG=false
SKIP_LOGS=false

if [[ "$*" == *"--noopt"* ]]; then
    RUN_NOOPT=true
fi

if [[ "$*" == *"--smart-prealigned"* ]]; then
    RUN_SMART_PREALIGNED=true
fi

if [[ "$*" == *"--approx-tail"* ]]; then
    RUN_APPROX_TAIL=true
fi

if [[ "$*" == *"--approx-dynamic"* ]]; then
    RUN_APPROX_DYNAMIC=true
fi

if [[ "$*" == *"--approx-adaptive-simple"* ]]; then
    RUN_APPROX_ADAPTIVE_SIMPLE=true
elif [[ "$*" == *"--approx-adaptive-gap"* ]]; then
    RUN_APPROX_ADAPTIVE_GAP=true
fi

if [[ "$*" == *"--approx-adaptive-single-shot"* ]]; then
    RUN_APPROX_ADAPTIVE_SINGLE_SHOT=true
fi



if [[ "$*" == *"--debug"* ]]; then
    PASH_DEBUG=true
fi

if [[ "$*" == *"--skip-logs"* ]]; then
    SKIP_LOGS=true
fi


NUM_REPEATS=1
if [[ "$*" =~ --repeats[[:space:]]+([0-9]+) ]]; then
    NUM_REPEATS="${BASH_REMATCH[1]}"
fi

# If no mode flags specified, run all modes (default behavior)
if [ "$RUN_NOOPT" = false ] && \
   [ "$RUN_SMART_PREALIGNED" = false ] && \
   [ "$RUN_APPROX_TAIL" = false ] && \
   [ "$RUN_APPROX_DYNAMIC" = false ] && \
   [ "$RUN_APPROX_ADAPTIVE_GAP" = false ] && \
   [ "$RUN_APPROX_ADAPTIVE_SIMPLE" = false ] && \
   [ "$RUN_APPROX_ADAPTIVE_SINGLE_SHOT" = false ]; then
    RUN_NOOPT=true
    RUN_SMART_PREALIGNED=true
    RUN_APPROX_TAIL=true
    RUN_APPROX_DYNAMIC=true
    RUN_APPROX_ADAPTIVE_GAP=true
    RUN_APPROX_ADAPTIVE_SIMPLE=true
    RUN_APPROX_ADAPTIVE_SINGLE_SHOT=true
fi

MODES=(
    noopt
    s3_smart_prealigned
    s3_approx_tail_coord
    s3_approx_dynamic
    s3_approx_adaptive_gap
    s3_approx_adaptive_simple
    s3_approx_adaptive_single_shot
)

# Reader strategy mapping:
#   s3_smart_prealigned   -> aws/s3-chunk-reader-smart-prealigned.py
#   s3_approx_tail_coord  -> aws/s3-chunk-reader-approx-tail-coordination.py
#   s3_approx_* (others)  -> aws/s3-chunk-reader-approx-correction.py
declare -A MODE_DESC MODE_ENV MODE_SUFFIX MODE_ENABLE_S3 MODE_ENABLED MODE_FLAG MODE_IS_BASELINE
declare -A MODE_TIMES MODE_BILLED_MS MODE_COST MODE_MATCH MODE_SPEEDUP MODE_COST_DIFF MODE_DIFF_EXCERPT MODE_LOCAL_FILE MODE_REP1_TIME
declare -A MODE_BILLED_MS_LIST MODE_COST_LIST

MODE_DESC[noopt]="WITHOUT S3 direct streaming optimization"
MODE_DESC[s3_smart_prealigned]="WITH S3 direct streaming - SMART prealigned chunks (EC2 boundary scan)"
MODE_DESC[s3_approx_tail_coord]="WITH S3 direct streaming - APPROX chunks + tail coordination (legacy)"
MODE_DESC[s3_approx_dynamic]="WITH S3 direct streaming - APPROX chunks + dynamic correction windows"
MODE_DESC[s3_approx_adaptive_gap]="WITH S3 direct streaming - APPROX chunks + adaptive gap-window (EC2-side)"
MODE_DESC[s3_approx_adaptive_simple]="WITH S3 direct streaming - APPROX chunks + adaptive simple (fixed sampled window)"
MODE_DESC[s3_approx_adaptive_single_shot]="WITH S3 direct streaming - APPROX chunks + adaptive single-shot sampled window"

MODE_FLAG[noopt]="--noopt"
MODE_FLAG[s3_smart_prealigned]="--smart-prealigned"
MODE_FLAG[s3_approx_tail_coord]="--approx-tail"
MODE_FLAG[s3_approx_dynamic]="--approx-dynamic"
MODE_FLAG[s3_approx_adaptive_gap]="--approx-adaptive-gap"
MODE_FLAG[s3_approx_adaptive_simple]="--approx-adaptive-simple"
MODE_FLAG[s3_approx_adaptive_single_shot]="--approx-adaptive-single-shot"

MODE_ENV[noopt]=""
MODE_ENV[s3_smart_prealigned]="USE_SMART_BOUNDARIES=true PASH_S3_CHUNKS_PER_LAMBDA=16"
MODE_ENV[s3_approx_tail_coord]="USE_SMART_BOUNDARIES=false"
MODE_ENV[s3_approx_dynamic]="PASH_S3_CHUNKS_PER_LAMBDA=16 USE_DYNAMIC_BOUNDARIES=true"
MODE_ENV[s3_approx_adaptive_gap]="USE_ADAPTIVE_BOUNDARIES=true PASH_GAP_SAMPLE_KB=256 PASH_GAP_DELTA=0.001 PASH_GAP_K_SAMPLES=4096 PASH_GAP_SAFETY_FACTOR=1.2 PASH_GAP_MAX_WINDOW_KB=1024 PASH_S3_CHUNKS_PER_LAMBDA=16"
MODE_ENV[s3_approx_adaptive_simple]="USE_ADAPTIVE_SIMPLE=true PASH_ADAPTIVE_SIMPLE_NUM_SAMPLES=5 PASH_ADAPTIVE_SIMPLE_SAMPLE_KB=256 PASH_ADAPTIVE_SIMPLE_SAFETY_FACTOR=1.5 PASH_S3_CHUNKS_PER_LAMBDA=16"
MODE_ENV[s3_approx_adaptive_single_shot]="USE_SINGLE_SHOT=true PASH_SINGLE_SHOT_SAMPLE_KB=256 PASH_SINGLE_SHOT_SAFETY_FACTOR=2.0 PASH_S3_CHUNKS_PER_LAMBDA=16"

MODE_SUFFIX[noopt]="noopt"
MODE_SUFFIX[s3_smart_prealigned]="s3smartprealigned"
MODE_SUFFIX[s3_approx_tail_coord]="s3approxtailcoord"
MODE_SUFFIX[s3_approx_dynamic]="s3approxdynamic"
MODE_SUFFIX[s3_approx_adaptive_gap]="s3approxadaptivegap"
MODE_SUFFIX[s3_approx_adaptive_simple]="s3approxadaptivesimple"
MODE_SUFFIX[s3_approx_adaptive_single_shot]="s3approxadaptivesingleshot"

MODE_ENABLE_S3[noopt]="false"
MODE_ENABLE_S3[s3_smart_prealigned]="true"
MODE_ENABLE_S3[s3_approx_tail_coord]="true"
MODE_ENABLE_S3[s3_approx_dynamic]="true"
MODE_ENABLE_S3[s3_approx_adaptive_gap]="true"
MODE_ENABLE_S3[s3_approx_adaptive_simple]="true"
MODE_ENABLE_S3[s3_approx_adaptive_single_shot]="true"

MODE_ENABLED[noopt]="$RUN_NOOPT"
MODE_ENABLED[s3_smart_prealigned]="$RUN_SMART_PREALIGNED"
MODE_ENABLED[s3_approx_tail_coord]="$RUN_APPROX_TAIL"
MODE_ENABLED[s3_approx_dynamic]="$RUN_APPROX_DYNAMIC"
MODE_ENABLED[s3_approx_adaptive_gap]="$RUN_APPROX_ADAPTIVE_GAP"
MODE_ENABLED[s3_approx_adaptive_simple]="$RUN_APPROX_ADAPTIVE_SIMPLE"
MODE_ENABLED[s3_approx_adaptive_single_shot]="$RUN_APPROX_ADAPTIVE_SINGLE_SHOT"

MODE_IS_BASELINE[noopt]="true"

# Benchmark-specific input matrix selection lives here.
source "${BENCHMARK_INPUT_CONFIG}"
if ! declare -F set_leash_benchmark_inputs >/dev/null; then
    echo "Error: ${BENCHMARK_INPUT_CONFIG} must define set_leash_benchmark_inputs()" >&2
    exit 2
fi

set_leash_benchmark_inputs "$@"
if [ "${#SCRIPT_INPUT_WIDTH[@]}" -eq 0 ]; then
    echo "Error: benchmark '${BENCHMARK_NAME}' produced empty SCRIPT_INPUT_WIDTH" >&2
    exit 2
fi

BENCHMARK_DIR="${BENCHMARK_NAME}"

# Check AWS_BUCKET is set
if [ -z "${AWS_BUCKET:-}" ]; then
    echo "Error: AWS_BUCKET environment variable is not set"
    exit 1
fi

if [ -z "${PASH_TOP:-}" ]; then
    echo "Error: PASH_TOP environment variable is not set"
    exit 1
fi

# Get lambda config for CSV output
read -r LAMBDA_MEM_MB LAMBDA_STORAGE_MB <<< "$(
  aws lambda get-function-configuration \
    --function-name lambda \
    --query '[MemorySize, EphemeralStorage.Size]' \
    --output text
)"

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

# Run pa.sh with timing and store LAST_WALL_TIME
run_pash_with_timing() {
    local mode_env="$1"
    local enable_s3="$2"
    local out_prefix="$3"
    local start_ns
    local end_ns

    start_ns=$(date +%s%N)
    if [ "$enable_s3" = "true" ]; then
        env PASH_DEBUG=$PASH_DEBUG $mode_env IN="$BENCHMARK_DIR/inputs/$INPUT" OUT="$out_prefix" DICT="oneliners/inputs/dict.txt" \
            $PASH_TOP/pa.sh --serverless_exec --enable_s3_direct -w"$WIDTH" scripts/"$SCRIPT"
    else
        env PASH_DEBUG=$PASH_DEBUG $mode_env IN="$BENCHMARK_DIR/inputs/$INPUT" OUT="$out_prefix" DICT="oneliners/inputs/dict.txt" \
            $PASH_TOP/pa.sh --serverless_exec -w"$WIDTH" scripts/"$SCRIPT"
    fi
    end_ns=$(date +%s%N)

    LAST_WALL_TIME=$(awk -v start="$start_ns" -v end="$end_ns" 'BEGIN{printf "%.3f", (end-start)/1000000000}')
}

# Fetch logs and parse billed duration/cost into LAST_BILLED_MS/LAST_COST
fetch_logs_and_cost() {
    local logs_dir="$1"
    local start_time_ms="$2"
    local job_id="$3"

    LAST_BILLED_MS="N/A"
    LAST_COST="N/A"

    if [ "$SKIP_LOGS" = false ]; then
        sleep 12

        if [ -d "$logs_dir" ]; then
            echo "Removing existing logs directory: $logs_dir"
            # du -sh "$logs_dir"
            # find "$logs_dir" -type f | wc -l
            # time rm -rf "$logs_dir"
            rm -rf "$logs_dir"
        fi
    
        local utils_output
        utils_output=$(python3 $PASH_TOP/scripts/serverless/utils.py "$logs_dir" "$start_time_ms" "$job_id")
        echo "$utils_output"

        LAST_BILLED_MS=$(echo "$utils_output" | sed -n 's/^\[Analysis\] Total billed time: \([0-9]\+\) ms/\1/p')
        LAST_COST=$(echo "$utils_output" | sed -n 's/^\[Analysis\] Cost estimate: \$\([0-9.]*\)/\1/p')
        if [ -z "$LAST_BILLED_MS" ]; then
            LAST_BILLED_MS="N/A"
        fi
        if [ -z "$LAST_COST" ]; then
            LAST_COST="N/A"
        fi
    else
        echo "Skipping log fetch (--skip-logs enabled)"
        sleep 1
    fi
}

# Download output and store LAST_LOCAL_FILE
download_mode_output() {
    local mode_suffix="$1"
    local out_prefix="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:${mode_suffix}"
    local s3_key="${out_prefix}stdout.txt"
    local local_file="/tmp/compare_${mode_suffix}_${SCRIPT//\//_}_${INPUT}_${WIDTH}.txt"

    LAST_LOCAL_FILE="$local_file"
    download_s3_output "$s3_key" "$local_file"
}

# Generic runner for modes (baseline included)
run_mode() {
    local mode="$1"
    local mode_index="$2"
    local mode_total="$3"
    local mode_suffix="${MODE_SUFFIX[$mode]}"
    local mode_desc="${MODE_DESC[$mode]}"
    local mode_env="${MODE_ENV[$mode]}"
    local enable_s3="${MODE_ENABLE_S3[$mode]}"
    local repeats="$NUM_REPEATS"
    local is_baseline="${MODE_IS_BASELINE[$mode]:-false}"

    if [ "$is_baseline" = "true" ]; then
        repeats=1
    fi

    echo "------------------------------------------------------------------------"
    echo "[MODE ${mode_index}/${mode_total}] ${mode_desc}"
    echo "Running $SCRIPT with input $INPUT and width $WIDTH"
    echo "------------------------------------------------------------------------"

    for REP in $(seq 1 "$repeats"); do
        echo "--- [$mode] Repeat $REP / $repeats ---"

        JOB_ID="${mode_suffix}_${SCRIPT//\//_}_${INPUT//\//_}_${WIDTH}_$(date +%s%N)_rep${REP}"
        export PASH_JOB_ID="$JOB_ID"
        echo "Job ID: $JOB_ID"

        START_TIME_MS=$(($(date +%s%3N) - 10000))
        echo "Start timestamp: $START_TIME_MS (with 10s safety buffer for clock skew)"

        local out_prefix="$BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:${mode_suffix}"
        local mode_wall_time

        run_pash_with_timing "$mode_env" "$enable_s3" "$out_prefix"
        mode_wall_time="$LAST_WALL_TIME"
        echo "[TIMING] ${mode} wall time: ${mode_wall_time}s"

        if [ -z "${MODE_TIMES[$mode]}" ]; then
            MODE_TIMES[$mode]="$mode_wall_time"
        else
            MODE_TIMES[$mode]="${MODE_TIMES[$mode]} $mode_wall_time"
        fi

        fetch_logs_and_cost "logs/$SCRIPT:$INPUT:$WIDTH:${mode_suffix}" "$START_TIME_MS" "$JOB_ID"
        local rep_billed_ms="$LAST_BILLED_MS"
        local rep_cost="$LAST_COST"

        if [ -z "${MODE_BILLED_MS_LIST[$mode]}" ]; then
            MODE_BILLED_MS_LIST[$mode]="$rep_billed_ms"
        else
            MODE_BILLED_MS_LIST[$mode]="${MODE_BILLED_MS_LIST[$mode]} $rep_billed_ms"
        fi
        if [ -z "${MODE_COST_LIST[$mode]}" ]; then
            MODE_COST_LIST[$mode]="$rep_cost"
        else
            MODE_COST_LIST[$mode]="${MODE_COST_LIST[$mode]} $rep_cost"
        fi

        if [ "$REP" -eq 1 ]; then
            MODE_REP1_TIME[$mode]="$mode_wall_time"
            MODE_BILLED_MS[$mode]="$rep_billed_ms"
            MODE_COST[$mode]="$rep_cost"

            if [ "$is_baseline" = "true" ]; then
                echo ""
                echo "Downloading ${mode_suffix} output from S3..."
                download_mode_output "$mode_suffix"
                MODE_LOCAL_FILE[$mode]="$LAST_LOCAL_FILE"
                echo ""
                NOOPT_WALL_TIME="$mode_wall_time"
                NOOPT_BILLED_MS="${MODE_BILLED_MS[$mode]}"
                NOOPT_COST="${MODE_COST[$mode]}"
                noopt_local_file="${MODE_LOCAL_FILE[$mode]}"
                write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "1" "$mode_wall_time" "$rep_billed_ms" "$rep_cost" "baseline" "baseline" "baseline" "baseline" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB"
            else
                MODE_SPEEDUP[$mode]="N/A"
                # Speedup calculation disabled per request.
                # if [[ "$NOOPT_WALL_TIME" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                #     MODE_SPEEDUP[$mode]=$(awk -v noopt="$NOOPT_WALL_TIME" -v mode="$mode_wall_time" 'BEGIN{ if (noopt>0 && mode>0) printf "%.3f", noopt/mode; else print "N/A" }')
                # fi

                MODE_COST_DIFF[$mode]="N/A"
                # Cost diff calculation disabled per request.
                # if [[ "${MODE_COST[$mode]}" =~ ^[0-9]+([.][0-9]+)?$ ]] && [[ "$NOOPT_COST" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                #     MODE_COST_DIFF[$mode]=$(awk -v mode="${MODE_COST[$mode]}" -v noopt="$NOOPT_COST" 'BEGIN{printf "%.6f", mode-noopt}')
                # fi
            fi
        else
            if [ "$is_baseline" = "true" ]; then
                echo "[REP $REP] Unexpected extra baseline repeat; skipping"
                continue
            fi

            echo "[REP $REP] Skipping output download"
            local rep_speedup="N/A"
            # Speedup calculation disabled per request.
            # if [[ "$NOOPT_WALL_TIME" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
            #     rep_speedup=$(awk -v noopt="$NOOPT_WALL_TIME" -v mode="$mode_wall_time" 'BEGIN{ if (noopt>0 && mode>0) printf "%.3f", noopt/mode; else print "N/A" }')
            # fi
            local rep_cost_diff="N/A"
            # Cost diff calculation disabled per request.
            # if [[ "$rep_cost" =~ ^[0-9]+([.][0-9]+)?$ ]] && [[ "$NOOPT_COST" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
            #     rep_cost_diff=$(awk -v mode="$rep_cost" -v noopt="$NOOPT_COST" 'BEGIN{printf "%.6f", mode-noopt}')
            # fi
            write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "$REP" "$mode_wall_time" "$rep_billed_ms" "$rep_cost" "$rep_speedup" "$rep_cost_diff" "N/A" "N/A" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB"
        fi
    done
}

# Helper function to compare two output files
compare_outputs() {
    local file1=$1
    local file2=$2
    local label1=$3
    local label2=$4
    DIFF_EXCERPT=""

    if [ ! -f "$file1" ]; then
        echo "✗ COMPARISON FAILED: $file1 does not exist"
        DIFF_EXCERPT="missing_file"
        return 1
    fi

    if [ ! -f "$file2" ]; then
        echo "✗ COMPARISON FAILED: $file2 does not exist"
        DIFF_EXCERPT="missing_file"
        return 1
    fi

    echo "Comparing outputs:"
    echo "  - $label1: $file1 ($(wc -l < "$file1") lines)"
    echo "  - $label2: $file2 ($(wc -l < "$file2") lines)"
    echo ""

    if cmp -s "$file1" "$file2" > /dev/null 2>&1; then
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "✓✓✓ OUTPUTS MATCH: Results are identical ✓✓✓"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        DIFF_EXCERPT=""
        return 0
    fi

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "✗✗✗ OUTPUTS DIFFER: Results are NOT identical ✗✗✗"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    return -1 # disable actual comparison for now since it will take long time for large files

    # If either file is empty, suppress unified diff headers/hunks
    local show_headers=1
    if [ ! -s "$file1" ] || [ ! -s "$file2" ]; then
        show_headers=0
    fi

    local tmpdiff
    tmpdiff="$(mktemp)"
    diff -u "$file1" "$file2" > "$tmpdiff" || true

    echo "Diff output (only changed lines; first 100 chars of each):"
    echo "────────────────────────────────────────────────────────────────────────"
    awk -v max=100 -v show_hdrs="$show_headers" '
        NR <= 2 { if (show_hdrs) print; next }
        /^@@/ { if (show_hdrs) print; next }
        /^[-+]/ {
            print substr($0, 1, max)
            next
        }

        # Skip context lines (those starting with space) and everything else
    ' "$tmpdiff"
    echo "────────────────────────────────────────────────────────────────────────"
    echo ""

    DIFF_EXCERPT="$(awk -v max=100 '
        NR <= 2 { next }
        /^@@/ { next }
        /^[-+]/ {
            line = substr($0, 1, max)
            gsub(/[;,|]/, " ", line)
            gsub(/\r/, "", line)
            if (out == "") out = line; else out = out ";" line
        }
        END { print out }
    ' "$tmpdiff")"

    local diff_lines
    diff_lines="$(awk 'NR > 2 && /^[-+]/ {c++} END{print c+0}' "$tmpdiff")"
    rm -f "$tmpdiff"

    echo "Total changed lines (additions+deletions): $diff_lines"
    return 1
}


# Array to track comparison results
comparison_results=()

# CSV helper and preamble
write_csv_row() {
    local run_start_time="$1"
    local script="$2"
    local input="$3"
    local width="$4"
    local mode="$5"
    local run_number="$6"
    local wall_time_sec="$7"
    local billed_duration_ms="$8"
    local cost_usd="$9"
    local benchmark="$BENCHMARK_NAME"
    # local speedup_vs_noopt="${10}"
    # local cost_diff_vs_noopt="${11}"
    local output_matches_noopt="${12}"
    local diff_excerpt="${13}"
    local lambda_mem_mb="${14}"
    local lambda_storage_mb="${15}"

    # echo "${run_start_time},${benchmark},${script},${input},${width},${mode},${run_number},${wall_time_sec},${billed_duration_ms},${cost_usd},${speedup_vs_noopt},${cost_diff_vs_noopt},${output_matches_noopt},${diff_excerpt}" >> "$RESULTS_CSV"
    echo "${run_start_time},${benchmark},${script},${input},${width},${mode},${run_number},${wall_time_sec},${billed_duration_ms},${cost_usd},${output_matches_noopt},${diff_excerpt},${lambda_mem_mb},${lambda_storage_mb}" >> "$RESULTS_CSV"
}

RUN_START_TIME=$(date +%Y-%m-%d_%H-%M-%S)
RESULTS_DIR="benchmark_results/$RUN_START_TIME"
mkdir -p "$RESULTS_DIR"
RESULTS_CSV="$RESULTS_DIR/results.csv"
# echo "run_start_time,benchmark,script,input,width,mode,run_number,wall_time_sec,billed_duration_ms,cost_usd,speedup_vs_noopt,cost_diff_vs_noopt,output_matches_noopt,diff_excerpt" > "$RESULTS_CSV"
echo "run_start_time,benchmark,script,input,width,mode,run_number,wall_time_sec,billed_duration_ms,cost_usd,output_matches_noopt,diff_excerpt,lambda_mem_mb,lambda_storage_mb" > "$RESULTS_CSV"
echo "CSV results: $RESULTS_CSV"

# Run benchmarks for all enabled modes
for SCRIPT_INPUT in "${SCRIPT_INPUT_WIDTH[@]}"; do
    echo "========================================================================"
    echo "Running benchmark for $SCRIPT_INPUT"
    SCRIPT=$(echo "$SCRIPT_INPUT" | cut -d: -f1)
    INPUT=$(echo "$SCRIPT_INPUT" | cut -d: -f2)
    WIDTH=$(echo "$SCRIPT_INPUT" | cut -d: -f3)
    WIDTH=${WIDTH:-64}

    NOOPT_WALL_TIME="N/A"
    NOOPT_BILLED_MS="N/A"
    NOOPT_COST="N/A"
    noopt_local_file=""
    for mode in "${MODES[@]}"; do
        MODE_TIMES[$mode]=""
        MODE_BILLED_MS[$mode]="N/A"
        MODE_COST[$mode]="N/A"
        MODE_MATCH[$mode]="N/A"
        MODE_SPEEDUP[$mode]="N/A"
        MODE_COST_DIFF[$mode]="N/A"
        MODE_DIFF_EXCERPT[$mode]=""
        MODE_LOCAL_FILE[$mode]=""
        MODE_REP1_TIME[$mode]="N/A"
        MODE_BILLED_MS_LIST[$mode]=""
        MODE_COST_LIST[$mode]=""
    done

    mode_total=${#MODES[@]}
    mode_index=1
    for mode in "${MODES[@]}"; do
        if [ "${MODE_ENABLED[$mode]}" = true ]; then
            run_mode "$mode" "$mode_index" "$mode_total"
        else
            echo "------------------------------------------------------------------------"
            echo "[MODE ${mode_index}/${mode_total}] SKIPPING ${mode} mode (${MODE_FLAG[$mode]} flag not specified)"
            echo "------------------------------------------------------------------------"
            echo ""
        fi
        mode_index=$((mode_index + 1))
    done

    # Compare outputs (baseline vs enabled modes)
    for mode in "${MODES[@]}"; do
        if [ "${MODE_IS_BASELINE[$mode]:-false}" = "true" ]; then
            continue
        fi
        if [ "${MODE_ENABLED[$mode]}" != true ]; then
            continue
        fi

        mode_suffix="${MODE_SUFFIX[$mode]}"

        echo ""
        echo "Downloading ${mode_suffix} output from S3..."
        download_mode_output "$mode_suffix"
        MODE_LOCAL_FILE[$mode]="$LAST_LOCAL_FILE"
        echo ""

        if [ "${MODE_ENABLED[noopt]}" = true ]; then
            echo ""
            echo "╔════════════════════════════════════════════════════════════════════════╗"
            echo "║  LIVE COMPARISON FOR: $SCRIPT_INPUT (noopt vs ${mode_suffix})"
            echo "╚════════════════════════════════════════════════════════════════════════╝"
            echo ""
            if compare_outputs "$noopt_local_file" "${MODE_LOCAL_FILE[$mode]}" "noopt" "${mode_suffix}"; then
                MODE_MATCH[$mode]="true"
                MODE_DIFF_EXCERPT[$mode]=""
                comparison_results+=("$SCRIPT_INPUT (${mode_suffix}): ✓ MATCH")
                echo ""
                echo "╔════════════════════════════════════════════════════════════════════════╗"
                echo "║  ✓ RESULT: OUTPUTS MATCH for $SCRIPT_INPUT (noopt vs ${mode_suffix})"
                echo "╚════════════════════════════════════════════════════════════════════════╝"
            else
                MODE_MATCH[$mode]="false"
                MODE_DIFF_EXCERPT[$mode]="$DIFF_EXCERPT"
                comparison_results+=("$SCRIPT_INPUT (${mode_suffix}): ✗ DIFFER")
                echo ""
                echo "╔════════════════════════════════════════════════════════════════════════╗"
                echo "║  ✗ RESULT: OUTPUTS DIFFER for $SCRIPT_INPUT (noopt vs ${mode_suffix})"
                echo "╚════════════════════════════════════════════════════════════════════════╝"
            fi

            write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "1" "${MODE_REP1_TIME[$mode]}" "${MODE_BILLED_MS[$mode]}" "${MODE_COST[$mode]}" "${MODE_SPEEDUP[$mode]}" "${MODE_COST_DIFF[$mode]}" "${MODE_MATCH[$mode]}" "${MODE_DIFF_EXCERPT[$mode]}" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB"
        else
            write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "1" "${MODE_REP1_TIME[$mode]}" "${MODE_BILLED_MS[$mode]}" "${MODE_COST[$mode]}" "N/A" "N/A" "N/A" "N/A" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB"
        fi

        echo "Removing ${mode_suffix} local file: ${MODE_LOCAL_FILE[$mode]}..."
        rm -f "${MODE_LOCAL_FILE[$mode]}"
    done

    echo "Removing noopt local file: $noopt_local_file..."
    rm -f "$noopt_local_file"

    if [ "$NUM_REPEATS" -gt 1 ]; then
        # Min/max/avg summary calculations disabled per request.
        # for mode in "${MODES[@]}"; do
        #     if [ "${MODE_IS_BASELINE[$mode]}" = "true" ]; then
        #         continue
        #     fi
        #     if [ "${MODE_ENABLED[$mode]}" != true ]; then
        #         continue
        #     fi
        #     if [ -z "${MODE_TIMES[$mode]}" ]; then
        #         continue
        #     fi
        #
        #     mode_suffix="${MODE_SUFFIX[$mode]}"
        #
        #     min_time=$(echo "${MODE_TIMES[$mode]}" | awk '{min=$1; for(i=1;i<=NF;i++) if($i<min) min=$i; printf "%.3f", min}')
        #     max_time=$(echo "${MODE_TIMES[$mode]}" | awk '{max=$1; for(i=1;i<=NF;i++) if($i>max) max=$i; printf "%.3f", max}')
        #     avg_time=$(echo "${MODE_TIMES[$mode]}" | awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; printf "%.3f", sum/NF}')
        #
        #     min_billed=$(echo "${MODE_BILLED_MS_LIST[$mode]}" | awk '{
        #         for(i=1;i<=NF;i++) if($i ~ /^[0-9]+([.][0-9]+)?$/) { if(min=="" || $i<min) min=$i }
        #     } END { if(min=="") print "N/A"; else printf "%.0f", min }')
        #     max_billed=$(echo "${MODE_BILLED_MS_LIST[$mode]}" | awk '{
        #         for(i=1;i<=NF;i++) if($i ~ /^[0-9]+([.][0-9]+)?$/) { if(max=="" || $i>max) max=$i }
        #     } END { if(max=="") print "N/A"; else printf "%.0f", max }')
        #     avg_billed=$(echo "${MODE_BILLED_MS_LIST[$mode]}" | awk '{
        #         for(i=1;i<=NF;i++) if($i ~ /^[0-9]+([.][0-9]+)?$/) { sum+=$i; count++ }
        #     } END { if(count==0) print "N/A"; else printf "%.0f", sum/count }')
        #
        #     min_cost=$(echo "${MODE_COST_LIST[$mode]}" | awk '{
        #         for(i=1;i<=NF;i++) if($i ~ /^[0-9]+([.][0-9]+)?$/) { if(min=="" || $i<min) min=$i }
        #     } END { if(min=="") print "N/A"; else printf "%.6f", min }')
        #     max_cost=$(echo "${MODE_COST_LIST[$mode]}" | awk '{
        #         for(i=1;i<=NF;i++) if($i ~ /^[0-9]+([.][0-9]+)?$/) { if(max=="" || $i>max) max=$i }
        #     } END { if(max=="") print "N/A"; else printf "%.6f", max }')
        #     avg_cost=$(echo "${MODE_COST_LIST[$mode]}" | awk '{
        #         for(i=1;i<=NF;i++) if($i ~ /^[0-9]+([.][0-9]+)?$/) { sum+=$i; count++ }
        #     } END { if(count==0) print "N/A"; else printf "%.6f", sum/count }')
        #
        #     avg_speedup="N/A"
        #     if [[ "$NOOPT_WALL_TIME" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        #         avg_speedup=$(awk -v noopt="$NOOPT_WALL_TIME" -v avg="$avg_time" 'BEGIN{ if (noopt>0 && avg>0) printf "%.3f", noopt/avg; else print "N/A" }')
        #     fi
        #
        #     avg_cost_diff="N/A"
        #     if [[ "$avg_cost" =~ ^[0-9]+([.][0-9]+)?$ ]] && [[ "$NOOPT_COST" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        #         avg_cost_diff=$(awk -v mode="$avg_cost" -v noopt="$NOOPT_COST" 'BEGIN{printf "%.6f", mode-noopt}')
        #     fi
        #
        #     write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "min" "$min_time" "$min_billed" "$min_cost" "N/A" "N/A" "N/A" "N/A"
        #     write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "max" "$max_time" "$max_billed" "$max_cost" "N/A" "N/A" "N/A" "N/A"
        #     write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "avg" "$avg_time" "$avg_billed" "$avg_cost" "$avg_speedup" "$avg_cost_diff" "N/A" "N/A"
        # done
        :
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
echo "  - Optimized modes:   $BENCHMARK_DIR/outputs/*:s3* and logs/*:s3*"
echo "  - Graphviz graphs:   pash_graphviz_noopt_* and pash_graphviz_s3_*"
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

echo ""
echo "CSV results: $RESULTS_CSV"
echo "CSV head:"
head -n 5 "$RESULTS_CSV"

# Clean up temporary comparison files
echo ""
echo "Cleaning up temporary files..."
rm -f /tmp/compare_*.txt
echo "Done!"
