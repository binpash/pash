#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$ROOT_DIR"

TMP_DIR="$(mktemp -d)"
cleanup() {
  if [[ -n "${READER_PID:-}" ]]; then
    kill "${READER_PID}" >/dev/null 2>&1 || true
  fi
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT

INPUT_FILE="${TMP_DIR}/input_stateful.txt"
EXPECTED_OUTPUT_FILE="${TMP_DIR}/expected_stateful.txt"
OUTPUT_FILE="${TMP_DIR}/output_stateful.txt"
RDV_KEY="ft-test-stateful-$(date +%s)-$RANDOM"

cat > "$INPUT_FILE" <<'EOF'
alpha-stateful
beta-stateful
gamma-stateful
EOF
cp "$INPUT_FILE" "$EXPECTED_OUTPUT_FILE"

# Receiver: read from TCP and write to stdout.
"${SCRIPT_DIR}/reader.sh" "$RDV_KEY" "recv-node" "send-node" > "$OUTPUT_FILE" &
READER_PID=$!
sleep 1

# Crash sender #1 after sending only partial payload.
{
  head -c 12 "$INPUT_FILE"
  tail -f /dev/null
} | "${SCRIPT_DIR}/sender.sh" "$RDV_KEY" "send-node" "recv-node" "false" 0 "job-stateful" "folder-stateful" "script-stateful" &
SENDER1_PID=$!
sleep 1
kill -9 "$SENDER1_PID" >/dev/null 2>&1 || true
wait "$SENDER1_PID" >/dev/null 2>&1 || true

# Crash sender #2 the same way.
{
  head -c 12 "$INPUT_FILE"
  tail -f /dev/null
} | "${SCRIPT_DIR}/sender.sh" "$RDV_KEY" "send-node" "recv-node" "false" 0 "job-stateful" "folder-stateful" "script-stateful" &
SENDER2_PID=$!
sleep 1
kill -9 "$SENDER2_PID" >/dev/null 2>&1 || true
wait "$SENDER2_PID" >/dev/null 2>&1 || true

# Final sender run: send full input and close stdin cleanly.
cat "$INPUT_FILE" | "${SCRIPT_DIR}/sender.sh" "$RDV_KEY" "send-node" "recv-node" "false" 0 "job-stateful" "folder-stateful" "script-stateful"

wait "$READER_PID"
READER_PID=""

if cmp -s "$EXPECTED_OUTPUT_FILE" "$OUTPUT_FILE"; then
  echo "test2_stateful: PASS"
else
  echo "test2_stateful: FAIL (output mismatch)" >&2
  exit 1
fi
