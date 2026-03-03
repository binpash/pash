#!/usr/bin/env bash
cd "$(dirname "$0")"
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <rdv_key> [me] [peer]" >&2
  exit 1
fi

RDV_KEY="$1"
ME="${2:-1}"
PEER="${3:-0}"

TMP_DIR="$(mktemp -d)"
FIFO_PATH="${TMP_DIR}/reader.out.fifo"
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

  rm -rf "$TMP_DIR"
  exit "$ec"
}
trap cleanup INT TERM EXIT

mkfifo "$FIFO_PATH"
ARG="recv*${RDV_KEY}*${ME}*${PEER}*${FIFO_PATH}"

echo "[reader.sh] Running pashlib-ft $ARG"

../target/release/pashlib-ft "$ARG" &
PASH_PID=$!

# Stream FIFO to stdout.
cat "$FIFO_PATH" &
CAT_PID=$!
wait "$CAT_PID"
wait "$PASH_PID"
