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
  --noopt, --noopt-no-resplitting, --smart-prealigned, --approx-tail,
  --approx-dynamic, --approx-adaptive-gap, --approx-adaptive-simple,
  --approx-adaptive-single-shot, --small/--medium/--large, --skip-logs, --debug, --repeats N, --no-hybrid, --hybrid,
  --verif / -v (skip running noopt; download existing noopt result from S3 for comparison; --noopt overrides)
  Short approx aliases: --approx-dyn, --approx-gap, --approx-simple, --approx-single, --approx-ss
  --time-to-crash N       Crash each lambda after N seconds (logs actual elapsed time)
  --chunk-start-idx N     Only process chunks with block_id >= N (default: 0)
  --lambda-crash-idx N    Apply --time-to-crash and --chunk-start-idx only to lambda index N (0-based)

Chunking flags (mutually exclusive; N can be comma-separated for sweeps, e.g. 1,2,4):
  --chunks-per-lambda N   Use a fixed N chunks per lambda (default: 16)
  --chunk-size [N]        Compute chunks-per-lambda from file size. N is in MB (default: 1)
                          Formula: ceil(F / (w * N * 1024 * 1024))

Width flag (comma-separated for sweeps, e.g. 32,64):
  --width N, -w N         Override width for all inputs (default: use width from leash-matrix-inputs.sh)

Design:
  - This script is the only runner implementation.
  - Each benchmark folder only defines SCRIPT_INPUT_WIDTH selection in:
      <benchmark>/leash-matrix-inputs.sh
EOF
            exit 0
            ;;
        --width|-w)
            RUNNER_ARGS+=("--width" "$2"); shift 2 ;;
        --width=*)
            RUNNER_ARGS+=("--width" "${1#*=}"); shift ;;
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
    echo "Error: benchmark folder is required (e.g. --benchmark unixfun)" >&2
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

benchmark_name_to_s3_input_name() {
    local name="$1"
    case "$name" in
        unixfun) echo "unix50" ;;
        covid) echo "covid-mts" ;;
        weather) echo "max-temp" ;;
        file-enc) echo "log-analysis" ;;
        analytics) echo "log-analysis" ;;
        *) echo "$name" ;;
    esac
}

benchmark_name_to_s3_output_name() {
    local name="$1"
    case "$name" in
        unixfun) echo "unix50" ;;
        covid) echo "covid-mts" ;;
        weather) echo "max-temp" ;;
        *) echo "$name" ;;
    esac
}

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

# Make Go-installed binaries available (needed for port-scan/zannotate)
export GOPATH="${GOPATH:-$HOME/go}"
export PATH="$PATH:$GOPATH/bin"
_go_local="${BENCHMARK_PATH}/go_install/go/bin"
[ -d "$_go_local" ] && export PATH="$PATH:$_go_local"
unset _go_local

normalize_runner_mode_aliases
set -- "${RUNNER_ARGS[@]}"

MODE_FLAGS=(
    --noopt
    --noopt-no-resplitting
    --smart-prealigned
    --approx-tail
    --approx-dynamic
    --approx-adaptive-gap
    --approx-adaptive-simple
    --approx-adaptive-single-shot
    --no-resplitting-approx-dynamic
    --noopt-no-hybrid
    --noopt-no-resplitting-no-hybrid
    --approx-dynamic-no-hybrid
    --no-resplitting-approx-dynamic-no-hybrid
)

ALL_ALLOWED_FLAGS=(
    "${MODE_FLAGS[@]}"
    --small
    --medium
    --large
    --skip-logs
    --debug
    --repeats
    --parallel_pipelines
    --parallel_pipelines_limit
    --chunks-per-lambda
    --chunk-size
    --width
    --time-to-crash
    --chunk-start-idx
    --lambda-crash-idx
    --no-hybrid
    --hybrid
    --verif
    -v
    --v
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
    local saw_chunks_per_lambda=false
    local saw_chunk_size=false

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

        if [[ "$arg" == "--time-to-crash" ]]; then
            if [ $((i + 1)) -ge "${#RUNNER_ARGS[@]}" ]; then
                echo "Error: --time-to-crash requires a positive integer value" >&2; exit 2
            fi
            local ttc_value="${RUNNER_ARGS[$((i + 1))]}"
            if ! [[ "$ttc_value" =~ ^[0-9]+$ ]] || [ "$ttc_value" -lt 1 ]; then
                echo "Error: --time-to-crash: '$ttc_value' must be a positive integer" >&2; exit 2
            fi
            i=$((i + 2)); continue
        fi

        if [[ "$arg" == "--chunk-start-idx" ]]; then
            if [ $((i + 1)) -ge "${#RUNNER_ARGS[@]}" ]; then
                echo "Error: --chunk-start-idx requires a non-negative integer value" >&2; exit 2
            fi
            local csi_value="${RUNNER_ARGS[$((i + 1))]}"
            if ! [[ "$csi_value" =~ ^[0-9]+$ ]]; then
                echo "Error: --chunk-start-idx: '$csi_value' must be a non-negative integer" >&2; exit 2
            fi
            i=$((i + 2)); continue
        fi

        if [[ "$arg" == "--lambda-crash-idx" ]]; then
            if [ $((i + 1)) -ge "${#RUNNER_ARGS[@]}" ]; then
                echo "Error: --lambda-crash-idx requires a non-negative integer value" >&2; exit 2
            fi
            local lci_value="${RUNNER_ARGS[$((i + 1))]}"
            if ! [[ "$lci_value" =~ ^[0-9]+$ ]]; then
                echo "Error: --lambda-crash-idx: '$lci_value' must be a non-negative integer" >&2; exit 2
            fi
            i=$((i + 2)); continue
        fi

        if [[ "$arg" == "--chunks-per-lambda" ]]; then
            saw_chunks_per_lambda=true
            if [ $((i + 1)) -ge "${#RUNNER_ARGS[@]}" ]; then
                echo "Error: --chunks-per-lambda requires a positive integer value" >&2
                exit 2
            fi
            local cpl_value="${RUNNER_ARGS[$((i + 1))]}"
            IFS=',' read -ra _cpl_tokens <<< "$cpl_value"
            for _tok in "${_cpl_tokens[@]}"; do
                if ! [[ "$_tok" =~ ^[0-9]+$ ]] || [ "$_tok" -lt 1 ]; then
                    echo "Error: --chunks-per-lambda: '$_tok' must be a positive integer" >&2; exit 2
                fi
            done
            i=$((i + 2))
            continue
        fi

        if [[ "$arg" == "--chunk-size" ]]; then
            saw_chunk_size=true
            if [ $((i + 1)) -lt "${#RUNNER_ARGS[@]}" ]; then
                local cs_next="${RUNNER_ARGS[$((i + 1))]}"
                if [[ "$cs_next" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
                    IFS=',' read -ra _cs_tokens <<< "$cs_next"
                    for _tok in "${_cs_tokens[@]}"; do
                        if [ "$_tok" -lt 1 ]; then
                            echo "Error: --chunk-size: '$_tok' must be >= 1" >&2; exit 2
                        fi
                    done
                    i=$((i + 2))
                    continue
                fi
            fi
            i=$((i + 1))
            continue
        fi

        if [[ "$arg" == "--width" ]]; then
            if [ $((i + 1)) -ge "${#RUNNER_ARGS[@]}" ]; then
                echo "Error: --width requires a value" >&2; exit 2
            fi
            local w_value="${RUNNER_ARGS[$((i + 1))]}"
            IFS=',' read -ra _w_tokens <<< "$w_value"
            for _tok in "${_w_tokens[@]}"; do
                if ! [[ "$_tok" =~ ^[0-9]+$ ]] || [ "$_tok" -lt 1 ]; then
                    echo "Error: --width: '$_tok' must be a positive integer" >&2; exit 2
                fi
            done
            i=$((i + 2)); continue
        fi

        if [[ "$arg" == --* ]] && ! is_allowed_runner_flag "$arg"; then
            echo "Error: mode not supported: $arg" >&2
            echo "Supported mode flags: ${MODE_FLAGS[*]}" >&2
            exit 2
        fi

        i=$((i + 1))
    done

    if [ "$saw_chunks_per_lambda" = true ] && [ "$saw_chunk_size" = true ]; then
        echo "Error: --chunks-per-lambda and --chunk-size are mutually exclusive" >&2
        exit 2
    fi
}

validate_runner_flags

# Usage:
#   --noopt                  : Baseline (EC2 split + Lambda compute)
#   --smart-prealigned       : Smart prealigned chunks (EC2 scans boundaries)
#   --approx-tail            : Approx chunks + Lambda tail coordination (legacy)
#   --approx-dynamic       : Approx chunks + dynamic correction window
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
RUN_NOOPT_NO_RESPLITTING=false
RUN_SMART_PREALIGNED=false
RUN_APPROX_TAIL=false
RUN_APPROX_DYNAMIC=false
RUN_APPROX_ADAPTIVE_GAP=false
RUN_APPROX_ADAPTIVE_SIMPLE=false
RUN_APPROX_ADAPTIVE_SINGLE_SHOT=false
PASH_DEBUG=false
SKIP_LOGS=false
PARALLEL_PIPELINES=false
PARALLEL_PIPELINES_LIMIT=""
RUN_APPROX_DYNAMIC_NO_RESPLITTING=false

if [[ " $* " == *" --noopt "* ]]; then
    RUN_NOOPT=true
fi

if [[ " $* " == *" --noopt-no-resplitting "* ]]; then
    RUN_NOOPT_NO_RESPLITTING=true
fi

if [[ "$*" == *"--smart-prealigned"* ]]; then
    RUN_SMART_PREALIGNED=true
fi

if [[ "$*" == *"--approx-tail"* ]]; then
    RUN_APPROX_TAIL=true
fi

if [[ " $* " == *" --approx-dynamic "* ]]; then
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

if [[ " $* " == *" --no-resplitting-approx-dynamic "* ]]; then
    RUN_APPROX_DYNAMIC_NO_RESPLITTING=true
fi

RUN_VERIF=false
if [[ "$*" == *"--verif"* ]] || [[ " $* " == *" -v "* ]]; then
    RUN_VERIF=true
fi

RUN_NO_HYBRID=false
if [[ " $* " == *" --no-hybrid "* ]]; then
    RUN_NO_HYBRID=true
fi

# Combined --*-no-hybrid convenience flags (imply --no-hybrid + the base mode flag)
if [[ " $* " == *" --noopt-no-hybrid "* ]]; then
    RUN_NOOPT=true
    RUN_NO_HYBRID=true
fi
if [[ " $* " == *" --noopt-no-resplitting-no-hybrid "* ]]; then
    RUN_NOOPT_NO_RESPLITTING=true
    RUN_NO_HYBRID=true
fi
if [[ " $* " == *" --approx-dynamic-no-hybrid "* ]]; then
    RUN_APPROX_DYNAMIC=true
    RUN_NO_HYBRID=true
fi
if [[ " $* " == *" --no-resplitting-approx-dynamic-no-hybrid "* ]]; then
    RUN_APPROX_DYNAMIC_NO_RESPLITTING=true
    RUN_NO_HYBRID=true
fi

RUN_HYBRID=false
if [[ " $* " == *" --hybrid "* ]]; then
    RUN_HYBRID=true
fi


if [[ "$*" == *"--parallel_pipelines"* ]]; then
    PARALLEL_PIPELINES=true
fi

if [[ "$*" == *"--parallel_pipelines_limit"* ]]; then
    if [[ "$*" =~ --parallel_pipelines_limit[[:space:]]+([0-9]+) ]]; then
        PARALLEL_PIPELINES_LIMIT="${BASH_REMATCH[1]}"
        echo "Parallel pipelines enabled with limit: $PARALLEL_PIPELINES_LIMIT"
    else
        echo "Error: --parallel_pipelines_limit requires a numeric value" >&2
        exit 2
    fi
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

TIME_TO_CRASH=""
if [[ "$*" =~ --time-to-crash[[:space:]]+([0-9]+) ]]; then
    TIME_TO_CRASH="${BASH_REMATCH[1]}"
fi

CHUNK_START_IDX=""
if [[ "$*" =~ --chunk-start-idx[[:space:]]+([0-9]+) ]]; then
    CHUNK_START_IDX="${BASH_REMATCH[1]}"
fi

LAMBDA_CRASH_IDX=""
if [[ "$*" =~ --lambda-crash-idx[[:space:]]+([0-9]+) ]]; then
    LAMBDA_CRASH_IDX="${BASH_REMATCH[1]}"
fi

CHUNKS_MODE="fixed"
CHUNKS_PER_LAMBDA_VALUES=(16)   # array; default single value
CHUNK_SIZE_MB_VALUES=(1)        # array; default single value
WIDTH_OVERRIDE_VALUES=()        # empty = use per-entry width from config
FIXED_CHUNKS_PER_LAMBDA=16      # set per sweep iteration
CHUNK_SIZE_MB=1                 # set per sweep iteration

if [[ "$*" =~ --chunks-per-lambda[[:space:]]+([^[:space:]]+) ]]; then
    IFS=',' read -ra CHUNKS_PER_LAMBDA_VALUES <<< "${BASH_REMATCH[1]}"
    CHUNKS_MODE="fixed"
fi

if [[ "$*" == *"--chunk-size"* ]]; then
    CHUNKS_MODE="dynamic"
    if [[ "$*" =~ --chunk-size[[:space:]]+([0-9][0-9,]*) ]]; then
        IFS=',' read -ra CHUNK_SIZE_MB_VALUES <<< "${BASH_REMATCH[1]}"
    fi
fi

if [[ "$*" =~ --width[[:space:]]+([^[:space:]]+) ]]; then
    IFS=',' read -ra WIDTH_OVERRIDE_VALUES <<< "${BASH_REMATCH[1]}"
fi

# If no mode flags specified, run all modes (default behavior)
if [ "$RUN_NOOPT" = false ] && \
   [ "$RUN_NOOPT_NO_RESPLITTING" = false ] && \
   [ "$RUN_SMART_PREALIGNED" = false ] && \
   [ "$RUN_APPROX_TAIL" = false ] && \
   [ "$RUN_APPROX_DYNAMIC" = false ] && \
   [ "$RUN_APPROX_ADAPTIVE_GAP" = false ] && \
   [ "$RUN_APPROX_ADAPTIVE_SIMPLE" = false ] && \
   [ "$RUN_APPROX_ADAPTIVE_SINGLE_SHOT" = false ] && \
   [ "$RUN_APPROX_DYNAMIC_NO_RESPLITTING" = false ]; then
    RUN_NOOPT=true
    RUN_NOOPT_NO_RESPLITTING=true
    RUN_SMART_PREALIGNED=true
    RUN_APPROX_TAIL=true
    RUN_APPROX_DYNAMIC=true
    RUN_APPROX_ADAPTIVE_GAP=true
    RUN_APPROX_ADAPTIVE_SIMPLE=true
    RUN_APPROX_ADAPTIVE_SINGLE_SHOT=true
fi

# --noopt overrides --verif: if both specified, noopt runs normally and verif is disabled
if [ "$RUN_VERIF" = "true" ] && [ "$RUN_NOOPT" = "true" ]; then
    RUN_VERIF=false
elif [ "$RUN_VERIF" = "true" ]; then
    RUN_NOOPT=false
fi

MODES=(
    noopt
    noopt_no_resplitting
    s3_smart_prealigned
    s3_approx_tail_coord
    s3_approx_dynamic
    s3_approx_adaptive_gap
    s3_approx_adaptive_simple
    s3_approx_adaptive_single_shot
    s3_approx_dynamic_no_resplitting
    noopt_no_hybrid
    noopt_no_resplitting_no_hybrid
    s3_approx_dynamic_no_hybrid
    s3_approx_dynamic_no_resplitting_no_hybrid
)

# Reader strategy mapping:
#   s3_smart_prealigned   -> aws/s3-chunk-reader-smart-prealigned.py
#   s3_approx_tail_coord  -> aws/s3-chunk-reader-approx-tail-coordination.py
#   s3_approx_* (others)  -> aws/s3-chunk-reader-approx-correction.py
declare -A MODE_DESC MODE_ENV MODE_SUFFIX MODE_ENABLE_S3 MODE_ENABLED MODE_FLAG MODE_IS_BASELINE MODE_USES_CHUNKS_PER_LAMBDA
declare -A MODE_TIMES MODE_BILLED_MS MODE_COST_LAMBDA MODE_COST_TOTAL MODE_MATCH MODE_SPEEDUP MODE_COST_DIFF MODE_DIFF_EXCERPT MODE_LOCAL_FILE MODE_REP1_TIME
declare -A MODE_BILLED_MS_LIST MODE_COST_LIST
declare -A EC2_PRICE

MODE_DESC[noopt]="WITHOUT S3 direct streaming optimization"
MODE_DESC[noopt_no_resplitting]="WITHOUT S3 direct streaming optimization, WITHOUT resplitting"
MODE_DESC[s3_smart_prealigned]="WITH S3 direct streaming - SMART prealigned chunks (EC2 boundary scan)"
MODE_DESC[s3_approx_tail_coord]="WITH S3 direct streaming - APPROX chunks + tail coordination (legacy)"
MODE_DESC[s3_approx_dynamic]="WITH S3 direct streaming - APPROX chunks + dynamic correction windows"
MODE_DESC[s3_approx_adaptive_gap]="WITH S3 direct streaming - APPROX chunks + adaptive gap-window (EC2-side)"
MODE_DESC[s3_approx_adaptive_simple]="WITH S3 direct streaming - APPROX chunks + adaptive simple (fixed sampled window)"
MODE_DESC[s3_approx_adaptive_single_shot]="WITH S3 direct streaming - APPROX chunks + adaptive single-shot sampled window"
MODE_DESC[s3_approx_dynamic_no_resplitting]="WITH S3 direct streaming - APPROX chunks + dynamic correction windows, WITHOUT resplitting and streaming to lambdas after direct s3"
MODE_DESC[noopt_no_hybrid]="WITHOUT S3 direct streaming optimization, no hybrid"
MODE_DESC[noopt_no_resplitting_no_hybrid]="WITHOUT S3 direct streaming optimization, WITHOUT resplitting, no hybrid"
MODE_DESC[s3_approx_dynamic_no_hybrid]="WITH S3 direct streaming - APPROX chunks + dynamic correction windows, no hybrid"
MODE_DESC[s3_approx_dynamic_no_resplitting_no_hybrid]="WITH S3 direct streaming - APPROX chunks + dynamic correction windows, WITHOUT resplitting and streaming to lambdas after direct s3, no hybrid"

MODE_FLAG[noopt]="--noopt"
MODE_FLAG[noopt_no_resplitting]="--noopt-no-resplitting --no_resplitting --ec2_width $(nproc)"
MODE_FLAG[s3_smart_prealigned]="--smart-prealigned"
MODE_FLAG[s3_approx_tail_coord]="--approx-tail"
MODE_FLAG[s3_approx_dynamic]="--approx-dynamic"
MODE_FLAG[s3_approx_adaptive_gap]="--approx-adaptive-gap"
MODE_FLAG[s3_approx_adaptive_simple]="--approx-adaptive-simple"
MODE_FLAG[s3_approx_adaptive_single_shot]="--approx-adaptive-single-shot"
MODE_FLAG[s3_approx_dynamic_no_resplitting]="--approx-dynamic --no_resplitting --ec2_width $(nproc)"
MODE_FLAG[noopt_no_hybrid]="--noopt"
MODE_FLAG[noopt_no_resplitting_no_hybrid]="--noopt-no-resplitting --no_resplitting --ec2_width $(nproc)"
MODE_FLAG[s3_approx_dynamic_no_hybrid]="--approx-dynamic"
MODE_FLAG[s3_approx_dynamic_no_resplitting_no_hybrid]="--approx-dynamic --no_resplitting --ec2_width $(nproc)"

MODE_ENV[noopt]="LEASH_DISABLE_PASHLIB_FT=true"
MODE_ENV[noopt_no_resplitting]="LEASH_DISABLE_PASHLIB_FT=true NO_RESPLITTING=true"
MODE_ENV[s3_smart_prealigned]="USE_SMART_BOUNDARIES=true"
MODE_ENV[s3_approx_tail_coord]="USE_SMART_BOUNDARIES=false"
MODE_ENV[s3_approx_dynamic]="USE_DYNAMIC_BOUNDARIES=true"
MODE_ENV[s3_approx_adaptive_gap]="USE_ADAPTIVE_BOUNDARIES=true PASH_GAP_SAMPLE_KB=256 PASH_GAP_DELTA=0.001 PASH_GAP_K_SAMPLES=4096 PASH_GAP_SAFETY_FACTOR=1.2 PASH_GAP_MAX_WINDOW_KB=1024"
MODE_ENV[s3_approx_adaptive_simple]="USE_ADAPTIVE_SIMPLE=true PASH_ADAPTIVE_SIMPLE_NUM_SAMPLES=5 PASH_ADAPTIVE_SIMPLE_SAMPLE_KB=256 PASH_ADAPTIVE_SIMPLE_SAFETY_FACTOR=1.5"
MODE_ENV[s3_approx_adaptive_single_shot]="USE_SINGLE_SHOT=true PASH_SINGLE_SHOT_SAMPLE_KB=256 PASH_SINGLE_SHOT_SAFETY_FACTOR=2.0"
MODE_ENV[s3_approx_dynamic_no_resplitting]="USE_DYNAMIC_BOUNDARIES=true NO_RESPLITTING=true"
MODE_ENV[noopt_no_hybrid]="LEASH_DISABLE_PASHLIB_FT=true"
MODE_ENV[noopt_no_resplitting_no_hybrid]="LEASH_DISABLE_PASHLIB_FT=true NO_RESPLITTING=true"
MODE_ENV[s3_approx_dynamic_no_hybrid]="USE_DYNAMIC_BOUNDARIES=true"
MODE_ENV[s3_approx_dynamic_no_resplitting_no_hybrid]="USE_DYNAMIC_BOUNDARIES=true NO_RESPLITTING=true"

MODE_USES_CHUNKS_PER_LAMBDA[noopt]="false"
MODE_USES_CHUNKS_PER_LAMBDA[noopt_no_resplitting]="false"
MODE_USES_CHUNKS_PER_LAMBDA[s3_smart_prealigned]="true"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_tail_coord]="false"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_dynamic]="true"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_adaptive_gap]="true"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_adaptive_simple]="true"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_adaptive_single_shot]="true"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_dynamic_no_resplitting]="true"
MODE_USES_CHUNKS_PER_LAMBDA[noopt_no_hybrid]="false"
MODE_USES_CHUNKS_PER_LAMBDA[noopt_no_resplitting_no_hybrid]="false"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_dynamic_no_hybrid]="true"
MODE_USES_CHUNKS_PER_LAMBDA[s3_approx_dynamic_no_resplitting_no_hybrid]="true"

MODE_SUFFIX[noopt]="noopt"
MODE_SUFFIX[noopt_no_resplitting]="nooptnoresplit"
MODE_SUFFIX[s3_smart_prealigned]="s3smartprealigned"
MODE_SUFFIX[s3_approx_tail_coord]="s3approxtailcoord"
MODE_SUFFIX[s3_approx_dynamic]="s3approxdynamic"
MODE_SUFFIX[s3_approx_adaptive_gap]="s3approxadaptivegap"
MODE_SUFFIX[s3_approx_adaptive_simple]="s3approxadaptivesimple"
MODE_SUFFIX[s3_approx_adaptive_single_shot]="s3approxadaptivesingleshot"
MODE_SUFFIX[s3_approx_dynamic_no_resplitting]="s3approxdynamicnoresplit"
MODE_SUFFIX[noopt_no_hybrid]="noopt_no_hybrid"
MODE_SUFFIX[noopt_no_resplitting_no_hybrid]="nooptnoresplit_no_hybrid"
MODE_SUFFIX[s3_approx_dynamic_no_hybrid]="s3approxdynamic_no_hybrid"
MODE_SUFFIX[s3_approx_dynamic_no_resplitting_no_hybrid]="s3approxdynamicnoresplit_no_hybrid"

MODE_ENABLE_S3[noopt]="false"
MODE_ENABLE_S3[noopt_no_resplitting]="false"
MODE_ENABLE_S3[s3_smart_prealigned]="true"
MODE_ENABLE_S3[s3_approx_tail_coord]="true"
MODE_ENABLE_S3[s3_approx_dynamic]="true"
MODE_ENABLE_S3[s3_approx_adaptive_gap]="true"
MODE_ENABLE_S3[s3_approx_adaptive_simple]="true"
MODE_ENABLE_S3[s3_approx_adaptive_single_shot]="true"
MODE_ENABLE_S3[s3_approx_dynamic_no_resplitting]="true"
MODE_ENABLE_S3[noopt_no_hybrid]="false"
MODE_ENABLE_S3[noopt_no_resplitting_no_hybrid]="false"
MODE_ENABLE_S3[s3_approx_dynamic_no_hybrid]="true"
MODE_ENABLE_S3[s3_approx_dynamic_no_resplitting_no_hybrid]="true"

MODE_ENABLED[noopt]="$RUN_NOOPT"
MODE_ENABLED[noopt_no_resplitting]="$RUN_NOOPT_NO_RESPLITTING"
MODE_ENABLED[s3_smart_prealigned]="$RUN_SMART_PREALIGNED"
MODE_ENABLED[s3_approx_tail_coord]="$RUN_APPROX_TAIL"
MODE_ENABLED[s3_approx_dynamic]="$RUN_APPROX_DYNAMIC"
MODE_ENABLED[s3_approx_adaptive_gap]="$RUN_APPROX_ADAPTIVE_GAP"
MODE_ENABLED[s3_approx_adaptive_simple]="$RUN_APPROX_ADAPTIVE_SIMPLE"
MODE_ENABLED[s3_approx_adaptive_single_shot]="$RUN_APPROX_ADAPTIVE_SINGLE_SHOT"
MODE_ENABLED[s3_approx_dynamic_no_resplitting]="$RUN_APPROX_DYNAMIC_NO_RESPLITTING"

# Adjust for hybrid/no-hybrid selection (only for the 4 no_hybrid-capable modes)
_NH_CAPABLE=(noopt noopt_no_resplitting s3_approx_dynamic s3_approx_dynamic_no_resplitting)
if [ "$RUN_NO_HYBRID" = "true" ] && [ "$RUN_HYBRID" != "true" ]; then
    # --no-hybrid only: enable no_hybrid variants, disable base variants
    for _nh_mode in "${_NH_CAPABLE[@]}"; do
        MODE_ENABLED[${_nh_mode}_no_hybrid]="${MODE_ENABLED[$_nh_mode]}"
        MODE_ENABLED[$_nh_mode]="false"
    done
elif [ "$RUN_NO_HYBRID" = "true" ] && [ "$RUN_HYBRID" = "true" ]; then
    # both --hybrid and --no-hybrid: keep base, also enable no_hybrid variants
    for _nh_mode in "${_NH_CAPABLE[@]}"; do
        MODE_ENABLED[${_nh_mode}_no_hybrid]="${MODE_ENABLED[$_nh_mode]}"
    done
else
    # neither --no-hybrid nor both: disable all _no_hybrid variants
    for _nh_mode in "${_NH_CAPABLE[@]}"; do
        MODE_ENABLED[${_nh_mode}_no_hybrid]="false"
    done
fi

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
S3_INPUT_BENCHMARK_DIR="$(benchmark_name_to_s3_input_name "${BENCHMARK_NAME}")"
S3_OUTPUT_BENCHMARK_DIR="$(benchmark_name_to_s3_output_name "${BENCHMARK_NAME}")"
echo "Benchmark: ${BENCHMARK_NAME}"
echo "Input benchmark S3 prefix: ${S3_INPUT_BENCHMARK_DIR}"
echo "Output benchmark S3 prefix: ${S3_OUTPUT_BENCHMARK_DIR}"

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

# Get ec2 config for CSV output and cost calculations
TOKEN=$(curl -s -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")
EC2_INSTANCE_TYPE=$(curl -s -H "X-aws-ec2-metadata-token: $TOKEN" "http://169.254.169.254/latest/meta-data/instance-type")

# Get EC2 base cost
EC2_PRICE["m5.large"]=0.096
EC2_PRICE["m5.xlarge"]=0.192
EC2_PRICE["m5.4xlarge"]=0.768
EC2_PRICE["m5.16large"]=3.072

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
    local no_resplitting_flag="$4"
    local no_hybrid_flag="$5"
    local entries_value="$6"
    local start_ns
    local end_ns

    start_ns=$(date +%s%N)
    parallel_config=""
    if [ "$PARALLEL_PIPELINES" = "true" ]; then
        parallel_config="--parallel_pipelines"
        if [ -n "$PARALLEL_PIPELINES_LIMIT" ]; then
            parallel_config+=" --parallel_pipelines_limit $PARALLEL_PIPELINES_LIMIT"
        fi
    fi
    benchmark_dir=$BENCHMARK_DIR
    if [[ $BENCHMARK_NAME == "file-enc" ]]; then
        benchmark_dir="analytics"
    fi
    if [ "$enable_s3" = "true" ]; then
        entries_value=${entries_value:-${LEASH_ENTRIES:-1}}
        env PASH_DEBUG=$PASH_DEBUG $mode_env IN="$S3_INPUT_BENCHMARK_DIR/inputs/$INPUT" OUT="$out_prefix" DICT="oneliners/inputs/dict.txt" ENTRIES=$entries_value \
            $PASH_TOP/pa.sh --serverless_exec --enable_s3_direct $no_resplitting_flag $no_hybrid_flag $parallel_config -w"$WIDTH" scripts/"$SCRIPT"
    else
        env PASH_DEBUG=$PASH_DEBUG $mode_env IN="$S3_INPUT_BENCHMARK_DIR/inputs/$INPUT" OUT="$out_prefix" DICT="oneliners/inputs/dict.txt" \
            $PASH_TOP/pa.sh --serverless_exec $no_resplitting_flag $no_hybrid_flag $parallel_config -w"$WIDTH" scripts/"$SCRIPT"
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
    local out_prefix="$S3_OUTPUT_BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:${mode_suffix}"

    if [ "$BENCHMARK_NAME" == "weather" ]; then
        local s3_key_1="${out_prefix}average.stdout.txt"
        local local_file_1="/tmp/compare_${mode_suffix}_${SCRIPT//\//_}_${INPUT}_${WIDTH}_average.txt"
        local s3_key_2="${out_prefix}min.stdout.txt"
        local local_file_2="/tmp/compare_${mode_suffix}_${SCRIPT//\//_}_${INPUT}_${WIDTH}_min.txt"
        local s3_key_3="${out_prefix}max.stdout.txt"
        local local_file_3="/tmp/compare_${mode_suffix}_${SCRIPT//\//_}_${INPUT}_${WIDTH}_max.txt"
        download_s3_output "$s3_key_1" "$local_file_1" || return 1
        download_s3_output "$s3_key_2" "$local_file_2" || return 1
        download_s3_output "$s3_key_3" "$local_file_3" || return 1
        local local_file="/tmp/compare_${mode_suffix}_${SCRIPT//\//_}_${INPUT}_${WIDTH}.txt"
        cat "$local_file_1" "$local_file_2" "$local_file_3" > "$local_file"
        rm "$local_file_1" "$local_file_2" "$local_file_3"
        LAST_LOCAL_FILE="$local_file"
        return 0
    fi

    local s3_key="${out_prefix}stdout.txt"
    local local_file="/tmp/compare_${mode_suffix}_${SCRIPT//\//_}_${INPUT}_${WIDTH}.txt"

    LAST_LOCAL_FILE="$local_file"
    download_s3_output "$s3_key" "$local_file"
}

# Benchmarks with large outputs where we intentionally skip file-by-file comparison.
should_skip_output_comparison() {
    [[ " nlp file-enc media-conv analytics web-search " == *" $BENCHMARK_NAME "* ]]
}

# Generic runner for modes (baseline included)
run_mode() {
    local mode="$1"
    local mode_index="$2"
    local mode_total="$3"
    local mode_suffix="${MODE_SUFFIX[$mode]}"
    local mode_desc="${MODE_DESC[$mode]}"
    local mode_env="${MODE_ENV[$mode]}"
    if [ "${MODE_USES_CHUNKS_PER_LAMBDA[$mode]}" = "true" ]; then
        mode_env="${mode_env:+$mode_env }PASH_S3_CHUNKS_PER_LAMBDA=${CURRENT_CHUNKS_PER_LAMBDA}"
    fi
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
        if [ -n "$TIME_TO_CRASH" ]; then
            export PASH_TIME_TO_CRASH="$TIME_TO_CRASH"
        fi
        if [ -n "$CHUNK_START_IDX" ]; then
            export PASH_CHUNK_START_IDX="$CHUNK_START_IDX"
        fi
        if [ -n "$LAMBDA_CRASH_IDX" ]; then
            export PASH_LAMBDA_CRASH_IDX="$LAMBDA_CRASH_IDX"
        fi

        START_TIME_MS=$(($(date +%s%3N) - 10000))
        echo "Start timestamp: $START_TIME_MS (with 10s safety buffer for clock skew)"

        local out_prefix="$S3_OUTPUT_BENCHMARK_DIR/outputs/$SCRIPT:$INPUT:$WIDTH:${mode_suffix}"
        local mode_wall_time

        no_resplitting=""
        if [[ "$mode" == "s3_approx_dynamic_no_resplitting" || "$mode" == "s3_approx_dynamic_no_resplitting_no_hybrid" ]]; then
            no_resplitting="--no_resplitting --ec2_width $(nproc)"
            if [[ " nlp file-enc media-conv analytics " == *" $BENCHMARK_NAME "* ]]; then
                no_resplitting="--no_resplitting --ec2_width 1 --unlimited_lambda"
            fi
            echo "Running APPROX DYNAMIC NO RESPLITTING mode $no_resplitting"
        elif [[ "$mode" == "noopt_no_resplitting" || "$mode" == "noopt_no_resplitting_no_hybrid" ]]; then
            no_resplitting="--no_resplitting --ec2_width $(nproc)"
        fi
        no_hybrid=""
        if [[ "$mode" == *"_no_hybrid" ]]; then
            no_hybrid="--no_hybrid"
        fi
        run_pash_with_timing "$mode_env" "$enable_s3" "$out_prefix" "$no_resplitting" "$no_hybrid" "$SCRIPT_LEASH_ENTRIES"
        mode_wall_time="$LAST_WALL_TIME"
        echo "[TIMING] ${mode} wall time: ${mode_wall_time}s"

        if [ -z "${MODE_TIMES[$mode]}" ]; then
            MODE_TIMES[$mode]="$mode_wall_time"
        else
            MODE_TIMES[$mode]="${MODE_TIMES[$mode]} $mode_wall_time"
        fi

        local ec2_cost="N/A"
        if [[ "${EC2_PRICE[$EC2_INSTANCE_TYPE]}" != "" ]] && [[ "$mode_wall_time" != "N/A" ]]; then
            ec2_cost=$(awk -v price="${EC2_PRICE[$EC2_INSTANCE_TYPE]}" -v time="$mode_wall_time" 'BEGIN{printf "%.6f", (price/3600)*time}')
        fi

        fetch_logs_and_cost "logs/$SCRIPT:$INPUT:$WIDTH:${mode_suffix}" "$START_TIME_MS" "$JOB_ID"
        local rep_billed_ms="$LAST_BILLED_MS"
        local rep_cost="$LAST_COST"

        local total_cost="N/A"
        if [[ "$ec2_cost" != "N/A" ]] && [[ "$LAST_COST" != "N/A" ]]; then
            total_cost=$(awk -v ec2="$ec2_cost" -v lambda="$LAST_COST" 'BEGIN{printf "%.6f", ec2+lambda}')
        fi

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
            MODE_COST_LAMBDA[$mode]="$rep_cost"
            MODE_COST_TOTAL[$mode]="$total_cost"

            if [ "$is_baseline" = "true" ]; then
                if ! should_skip_output_comparison; then
                    echo ""
                    echo "Downloading ${mode_suffix} output from S3..."
                    download_mode_output "$mode_suffix" || true
                    MODE_LOCAL_FILE[$mode]="${LAST_LOCAL_FILE:-}"
                    echo ""
                fi
                NOOPT_WALL_TIME="$mode_wall_time"
                NOOPT_BILLED_MS="${MODE_BILLED_MS[$mode]}"
                NOOPT_COST="${MODE_COST_LAMBDA[$mode]}"
                NOOPT_TOTAL_COST="${MODE_COST_TOTAL[$mode]}"
                noopt_local_file="${MODE_LOCAL_FILE[$mode]}"
                write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "1" "$mode_wall_time" "$rep_billed_ms" "$rep_cost" "$total_cost" "baseline" "baseline" "baseline" "baseline" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB" "$CURRENT_CHUNKS_PER_LAMBDA" "$CHUNK_SIZE_MB"
            else
                MODE_SPEEDUP[$mode]="N/A"
                # Speedup calculation disabled per request.
                # if [[ "$NOOPT_WALL_TIME" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                #     MODE_SPEEDUP[$mode]=$(awk -v noopt="$NOOPT_WALL_TIME" -v mode="$mode_wall_time" 'BEGIN{ if (noopt>0 && mode>0) printf "%.3f", noopt/mode; else print "N/A" }')
                # fi

                MODE_COST_DIFF[$mode]="N/A"
                # Cost diff calculation disabled per request.
                # if [[ "${MODE_COST_LAMBDA[$mode]}" =~ ^[0-9]+([.][0-9]+)?$ ]] && [[ "$NOOPT_COST" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                #     MODE_COST_DIFF[$mode]=$(awk -v mode="${MODE_COST_LAMBDA[$mode]}" -v noopt="$NOOPT_COST" 'BEGIN{printf "%.6f", mode-noopt}')
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
            write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "$REP" "$mode_wall_time" "$rep_billed_ms" "$rep_cost" "$total_cost" "$rep_speedup" "$rep_cost_diff" "N/A" "N/A" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB" "$CURRENT_CHUNKS_PER_LAMBDA" "$CHUNK_SIZE_MB"
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
    local cost_usd_lambda="$9"
    local cost_usd_total="${10}"
    local benchmark="$BENCHMARK_NAME"
    # local speedup_vs_noopt="${11}"
    # local cost_diff_vs_noopt="${12}"
    local output_matches_noopt="${13}"
    local diff_excerpt="${14}"
    local lambda_mem_mb="${15}"
    local lambda_storage_mb="${16}"
    local chunks_per_lambda="${17}"
    local chunk_size_mb="${18}"

    # echo "${run_start_time},${benchmark},${script},${input},${width},${mode},${run_number},${wall_time_sec},${billed_duration_ms},${cost_usd},${speedup_vs_noopt},${cost_diff_vs_noopt},${output_matches_noopt},${diff_excerpt}" >> "$RESULTS_CSV"
    echo "${run_start_time},${benchmark},${script},${input},${width},${mode},${run_number},${wall_time_sec},${billed_duration_ms},${cost_usd_lambda},${cost_usd_total},${output_matches_noopt},${diff_excerpt},${lambda_mem_mb},${lambda_storage_mb},${EC2_INSTANCE_TYPE},${chunks_per_lambda},${chunk_size_mb}" >> "$RESULTS_CSV"
}

RUN_START_TIME=$(date +%Y-%m-%d_%H-%M-%S)
RESULTS_DIR="benchmark_results/$RUN_START_TIME"
mkdir -p "$RESULTS_DIR"
RESULTS_CSV="$RESULTS_DIR/results.csv"
# echo "run_start_time,benchmark,script,input,width,mode,run_number,wall_time_sec,billed_duration_ms,cost_usd,speedup_vs_noopt,cost_diff_vs_noopt,output_matches_noopt,diff_excerpt" > "$RESULTS_CSV"
echo "run_start_time,benchmark,script,input,width,mode,run_number,wall_time_sec,billed_duration_ms,cost_usd_lambda,cost_usd_total,output_matches_noopt,diff_excerpt,lambda_mem_mb,lambda_storage_mb,ec2_config,chunks_per_lambda,chunk_size_mb" > "$RESULTS_CSV"
echo "CSV results: $RESULTS_CSV"

# Run benchmarks for all enabled modes

# Width sweep: if --width not given, use sentinel "" (per-entry width from config)
_WIDTH_SWEEP=("${WIDTH_OVERRIDE_VALUES[@]}")
[ "${#_WIDTH_SWEEP[@]}" -eq 0 ] && _WIDTH_SWEEP=("")

# Chunking sweep
if [ "$CHUNKS_MODE" = "fixed" ]; then
    _CHUNK_SWEEP=("${CHUNKS_PER_LAMBDA_VALUES[@]}")
else
    _CHUNK_SWEEP=("${CHUNK_SIZE_MB_VALUES[@]}")
fi

for _W in "${_WIDTH_SWEEP[@]}"; do
    [ -n "$_W" ] && echo "######## width sweep: ${_W} ########"

    for _C in "${_CHUNK_SWEEP[@]}"; do
        if [ "$CHUNKS_MODE" = "fixed" ]; then
            FIXED_CHUNKS_PER_LAMBDA="$_C"
            echo "######## chunks-per-lambda sweep: ${FIXED_CHUNKS_PER_LAMBDA} ########"
        else
            CHUNK_SIZE_MB="$_C"
            echo "######## chunk-size sweep: ${CHUNK_SIZE_MB} MB ########"
        fi

        for SCRIPT_INPUT in "${SCRIPT_INPUT_WIDTH[@]}"; do
    echo "========================================================================"
    echo "Running benchmark for $SCRIPT_INPUT"
    IFS=':' read -r SCRIPT INPUT WIDTH SCRIPT_LEASH_ENTRIES <<< "$SCRIPT_INPUT"
    WIDTH=${WIDTH:-64}
    [ -n "$_W" ] && WIDTH="$_W"

    _cs_file_size=$(aws s3api head-object \
        --bucket "$AWS_BUCKET" \
        --key "$S3_INPUT_BENCHMARK_DIR/inputs/$INPUT" \
        --query ContentLength --output text 2>/dev/null) || true

    if [ "$CHUNKS_MODE" = "fixed" ]; then
        CURRENT_CHUNKS_PER_LAMBDA="$FIXED_CHUNKS_PER_LAMBDA"
        if [[ "$_cs_file_size" =~ ^[0-9]+$ ]]; then
            CHUNK_SIZE_MB=$(awk "BEGIN{printf \"%.4g\", $_cs_file_size / ($WIDTH * $CURRENT_CHUNKS_PER_LAMBDA * 1024 * 1024)}")
        fi
    else
        if [[ "$_cs_file_size" =~ ^[0-9]+$ ]]; then
            _cs_chunk_bytes=$(awk "BEGIN{printf \"%d\", $CHUNK_SIZE_MB * 1024 * 1024}")
            _cs_denom=$(( WIDTH * _cs_chunk_bytes ))
            CURRENT_CHUNKS_PER_LAMBDA=$(( (_cs_file_size + _cs_denom - 1) / _cs_denom ))
            [ "$CURRENT_CHUNKS_PER_LAMBDA" -lt 1 ] && CURRENT_CHUNKS_PER_LAMBDA=1
        else
            echo "Warning: could not determine file size for '$INPUT'; defaulting chunks_per_lambda=16" >&2
            CURRENT_CHUNKS_PER_LAMBDA=16
        fi
    fi
    _cs_display=$(printf "%.4g" "$CHUNK_SIZE_MB")
    echo "[chunks] INPUT=$INPUT WIDTH=$WIDTH chunk_size_mb=${_cs_display} chunks_per_lambda=${CURRENT_CHUNKS_PER_LAMBDA}"

    NOOPT_WALL_TIME="N/A"
    NOOPT_BILLED_MS="N/A"
    NOOPT_COST="N/A"
    NOOPT_TOTAL_COST="N/A"
    noopt_local_file=""
    for mode in "${MODES[@]}"; do
        MODE_TIMES[$mode]=""
        MODE_BILLED_MS[$mode]="N/A"
        MODE_COST_LAMBDA[$mode]="N/A"
        MODE_COST_TOTAL[$mode]="N/A"
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

        mode_suffix="${MODE_SUFFIX[$mode]}"
        if [ "${MODE_ENABLED[$mode]}" != true ]; then
            continue
        fi

        if should_skip_output_comparison; then
            echo "Skipping output comparison for $BENCHMARK_NAME benchmark"
            write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "1" "${MODE_REP1_TIME[$mode]}" "${MODE_BILLED_MS[$mode]}" "${MODE_COST_LAMBDA[$mode]}" "${MODE_COST_TOTAL[$mode]}" "N/A" "N/A" "N/A" "N/A" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB" "$CURRENT_CHUNKS_PER_LAMBDA" "$CHUNK_SIZE_MB"
            continue
        fi

        echo ""
        echo "Downloading ${mode_suffix} output from S3..."
        download_mode_output "$mode_suffix" || true
        MODE_LOCAL_FILE[$mode]="${LAST_LOCAL_FILE:-}"
        echo ""

        # --verif: scan S3 for an existing noopt result matching same script+input (any width)
        if [ "$RUN_VERIF" = "true" ] && { [ -z "$noopt_local_file" ] || [ ! -f "$noopt_local_file" ]; }; then
            echo ""
            echo "Looking for existing noopt output on S3 (--verif mode, any width)..."
            noopt_s3_obj=$(aws s3 ls "s3://$AWS_BUCKET/$S3_OUTPUT_BENCHMARK_DIR/outputs/${SCRIPT}:${INPUT}:" 2>/dev/null \
                | awk '{print $NF}' \
                | grep -E "^${SCRIPT}:${INPUT}:[^:]+:nooptstdout\\.txt$" \
                | head -1)
            if [ -n "$noopt_s3_obj" ]; then
                echo "Found noopt result: $noopt_s3_obj"
                noopt_s3_key="$S3_OUTPUT_BENCHMARK_DIR/outputs/$noopt_s3_obj"
                noopt_verif_local="/tmp/compare_noopt_verif_${SCRIPT//\//_}_${INPUT}.txt"
                if download_s3_output "$noopt_s3_key" "$noopt_verif_local"; then
                    noopt_local_file="$noopt_verif_local"
                else
                    echo "✗ Could not download noopt output for verif comparison"
                fi
            else
                echo "✗ No noopt result found on S3 for $SCRIPT:$INPUT (any width)"
            fi
            echo ""
        fi

        if [ "${MODE_ENABLED[noopt]}" = true ] || [ "${MODE_ENABLED[noopt_no_hybrid]}" = true ] || [ -n "$noopt_local_file" ]; then
            if [ -z "$noopt_local_file" ] || [ ! -f "$noopt_local_file" ]; then
                echo ""
                echo "Downloading noopt output from S3..."
                download_mode_output "${MODE_SUFFIX[noopt]}" || true
                MODE_LOCAL_FILE[noopt]="${LAST_LOCAL_FILE:-}"
                noopt_local_file="${LAST_LOCAL_FILE:-}"
                echo ""
            fi

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

            write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "1" "${MODE_REP1_TIME[$mode]}" "${MODE_BILLED_MS[$mode]}" "${MODE_COST_LAMBDA[$mode]}" "${MODE_COST_TOTAL[$mode]}" "${MODE_SPEEDUP[$mode]}" "${MODE_COST_DIFF[$mode]}" "${MODE_MATCH[$mode]}" "${MODE_DIFF_EXCERPT[$mode]}" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB" "$CURRENT_CHUNKS_PER_LAMBDA" "$CHUNK_SIZE_MB"
        else
            write_csv_row "$RUN_START_TIME" "$SCRIPT" "$INPUT" "$WIDTH" "$mode_suffix" "1" "${MODE_REP1_TIME[$mode]}" "${MODE_BILLED_MS[$mode]}" "${MODE_COST_LAMBDA[$mode]}" "${MODE_COST_TOTAL[$mode]}" "N/A" "N/A" "N/A" "N/A" "$LAMBDA_MEM_MB" "$LAMBDA_STORAGE_MB" "$CURRENT_CHUNKS_PER_LAMBDA" "$CHUNK_SIZE_MB"
        fi

        echo "Removing ${mode_suffix} local file: ${MODE_LOCAL_FILE[$mode]}..."
        rm -f "${MODE_LOCAL_FILE[$mode]}"
    done

    if [ -n "$noopt_local_file" ]; then
        echo "Removing noopt local file: $noopt_local_file..."
        rm -f "$noopt_local_file"
    fi

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
        done  # end for SCRIPT_INPUT
    done  # end for _C (chunking sweep)
done  # end for _W (width sweep)

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
        match_count=$((match_count + 1))
        echo "  ✓ $result"
    else
        differ_count=$((differ_count + 1))
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
# rm -f /tmp/compare_*.txt
echo "Done!"
