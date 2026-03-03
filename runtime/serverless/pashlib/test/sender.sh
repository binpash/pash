#!/usr/bin/env bash
cd "$(dirname "$0")"
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <rdv_key> [me] [peer] [is_stateless] [chunk_start_id] [job_id] [folders_id] [script_id]" >&2
  exit 1
fi

RDV_KEY="$1"
ME="${2:-0}"
PEER="${3:-1}"
IS_STATELESS="${4:-false}"
CHUNK_START_ID="${5:-0}"
JOB_ID="${6:-test-job}"
FOLDERS_ID="${7:-test-folder}"
SCRIPT_ID="${8:-test-script}"

export LEASH_JOB_ID="$JOB_ID"
export FOLDERS_ID="$FOLDERS_ID"
export SCRIPT_ID="$SCRIPT_ID"
export CHUNK_START_ID="$CHUNK_START_ID"
export IS_STATELESS="$IS_STATELESS"

# make a random fifo
TEMP_DIR="$(mktemp -d)"
FIFO_PATH="${TEMP_DIR}/sender.in.fifo"
PASH_PID=""
CAT_PID=""

cleanup() {
  local ec=$?
  trap - INT TERM EXIT

  if [[ -n "${CAT_PID}" ]]; then
    kill -TERM "${CAT_PID}" >/dev/null 2>&1 || true
  fi
  if [[ -n "${PASH_PID}" ]]; then
    kill -TERM "${PASH_PID}" >/dev/null 2>&1 || true
    sleep 0.2
    kill -KILL "${PASH_PID}" >/dev/null 2>&1 || true
  fi

  rm -rf "$TEMP_DIR"
  exit "$ec"
}
trap cleanup INT TERM EXIT

mkfifo "$FIFO_PATH"
ARG="send*${RDV_KEY}*${ME}*${PEER}*${FIFO_PATH}"

echo "[sender.sh] Running pashlib-ft $ARG, is_stateless=$IS_STATELESS, chunk_start_id=$CHUNK_START_ID, job_id=$JOB_ID, folders_id=$FOLDERS_ID, script_id=$SCRIPT_ID"


# Stream stdin into FIFO; pashlib-ft reads from FIFO path.
cat > "$FIFO_PATH" &
CAT_PID=$!

../target/release/pashlib-ft "$ARG" &
PASH_PID=$!
wait "$PASH_PID"
wait "$CAT_PID"
