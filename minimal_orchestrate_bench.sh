#!/bin/bash
# Parameterized orchestration for holepunch benchmarking
# Usage: ./minimal_orchestrate_bench.sh <file_size> <output_mode>
# Example: ./minimal_orchestrate_bench.sh 100M /dev/null

set -e

export PASH_TOP="/home/ubuntu/pash"
INVOKE_SCRIPT="$PASH_TOP/aws/invoke-lambda.py"

# Parameters
FILE_SIZE=${1:-100M}        # 100M, 500M, or 1G
OUTPUT_MODE=${2:-/dev/null}  # /dev/null or /tmp/f.out

# Use timestamp as folder ID for unique S3 paths
FOLDER_ID=$(date +%s)

# Script identifier for Lambda
LAMBDA_SCRIPT_ID="lambda_bench"

# Start total timing
total_start=$(date +%s%3N)

# Delete old RDV entry (connection ID 42) to avoid stale data
aws dynamodb delete-item \
    --table-name rdv \
    --key '{"key":{"S":"42"}}' \
    --region us-east-1 >/dev/null 2>&1 || true

# Create Lambda script with inlined content (no external dependencies)
cat > /tmp/lambda_bench_${FOLDER_ID}.sh <<'EOFSCRIPT'
#!/bin/bash
export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1

# Benchmark parameters (substituted by orchestrator)
OUTPUT_MODE="OUTPUT_MODE_PLACEHOLDER"
version=$2

# Start timing (milliseconds since epoch)
start_time=$(date +%s%3N)

mkdir -p /tmp/pash_p4qiBDD/
mkdir -p /tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/
inform_daemon_exit ()
{
    msg="Exit:${process_id}";
    daemon_response=$(pash_communicate_daemon_just_send "$msg")
}
pash_communicate_daemon ()
{
    local message=$1;
    pash_communicate_unix_socket "compilation-server" "${DAEMON_SOCKET}" "${message}"
}
declare -fx pash_communicate_daemon
pash_communicate_daemon_just_send ()
{
    pash_communicate_daemon "$1"
}
declare -fx pash_communicate_daemon_just_send
pash_communicate_unix_socket ()
{
    local server_name=$1;
    local socket=$2;
    local message=$3;
    pash_redir_output echo "Sending msg to ${server_name}: $message";
    daemon_response=$(echo "$message" | nc -U "${socket}");
    pash_redir_output echo "Got response from ${server_name}: $daemon_response";
    echo "$daemon_response"
}
declare -fx pash_communicate_unix_socket
pash_redir_all_output ()
{
    :
}
declare -fx pash_redir_all_output
pash_redir_all_output_always_execute ()
{
    "$@" > /dev/null 2>&1
}
declare -fx pash_redir_all_output_always_execute
pash_redir_output ()
{
    :
}
declare -fx pash_redir_output
pash_wait_until_daemon_listening ()
{
    pash_wait_until_unix_socket_listening "compilation-server" "${DAEMON_SOCKET}"
}
declare -fx pash_wait_until_daemon_listening
pash_wait_until_unix_socket_listening ()
{
    local server_name=$1;
    local socket=$2;
    i=0;
    maximum_retries=1000;
    until echo "Daemon Start" 2> /dev/null | nc -U "$socket" > /dev/null 2>&1; do
        sleep 0.01;
        i=$((i+1));
        if [ $i -eq $maximum_retries ]; then
            echo "Error: Maximum retries: $maximum_retries exceeded when waiting for server: ${server_name} to bind to socket: ${socket}!" 1>&2;
            echo "Exiting..." 1>&2;
            exit 1;
        fi;
    done
}
declare -fx pash_wait_until_unix_socket_listening
run_parallel ()
{
    trap inform_daemon_exit SIGTERM SIGINT EXIT;
    export SCRIPT_TO_EXECUTE="$pash_script_to_execute";
    source "$RUNTIME_DIR/pash_restore_state_and_execute.sh"
}

rm_pash_fifos() {
{ rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo12" ; }
 { { rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo13" ; }
 { rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo18" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo12" ; }
 { { mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo13" ; }
 { mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo18" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""

# Consumer: use dgsh-tee for buffered I/O with parameterized output
{ runtime/dgsh-tee -i "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo12" -o "${OUTPUT_MODE}" -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"

{ /opt/pashlib recv*42*1*0*/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo12 & }
pids_to_kill="${!} ${pids_to_kill}"

source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}

# End timing and output
end_time=$(date +%s%3N)
duration_ms=$((end_time - start_time))
duration_sec=$(awk "BEGIN {printf \"%.3f\", $duration_ms/1000}")
echo "LAMBDA_TIME: ${duration_sec}s" >&2

rm_pash_fifos
( exit "${internal_exec_status}" )
EOFSCRIPT

# Substitute the OUTPUT_MODE placeholder with actual value
sed -i "s|OUTPUT_MODE_PLACEHOLDER|${OUTPUT_MODE}|g" /tmp/lambda_bench_${FOLDER_ID}.sh

aws s3 cp /tmp/lambda_bench_${FOLDER_ID}.sh \
    "s3://$AWS_BUCKET/sls-scripts/$FOLDER_ID/${LAMBDA_SCRIPT_ID}.sh" >/dev/null 2>&1

# Start both nodes in parallel
(cd "$PASH_TOP" && /bin/bash "$PASH_TOP/minimal/ec2_bench.sh" "$FILE_SIZE" "$FOLDER_ID" 2>&1) &
EC2_PID=$!

python3 "$INVOKE_SCRIPT" "$LAMBDA_SCRIPT_ID" "$FOLDER_ID" "lambda" >/dev/null 2>&1 &
LAMBDA_PID=$!

# Wait for both to complete
wait $EC2_PID $LAMBDA_PID
EXIT_CODE=$?

# Calculate total time
total_end=$(date +%s%3N)
total_duration_ms=$((total_end - total_start))
total_duration_sec=$(awk "BEGIN {printf \"%.3f\", $total_duration_ms/1000}")
echo "TOTAL_TIME: ${total_duration_sec}s" >&2

# Cleanup
rm -f /tmp/lambda_bench_${FOLDER_ID}.sh

exit $EXIT_CODE
