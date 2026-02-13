PARALLEL_PIPELINES_LIMITS=""
if [[ "$*" == *"--parallel_pipelines_limits"* ]]; then
    if [[ "$*" =~ --parallel_pipelines_limits[[:space:]]+([0-9]+) ]]; then
        PARALLEL_PIPELINES_LIMITS="${BASH_REMATCH[1]}"
        echo "Parallel pipelines enabled with limits: $PARALLEL_PIPELINES_LIMITS"
    else
        echo "Error: --parallel_pipelines_limits requires a numeric value" >&2
        exit 2
    fi
fi
if [ -n "$PARALLEL_PIPELINES_LIMITS" ]; then
    echo "$PARALLEL_PIPELINES_LIMITS"
fi
