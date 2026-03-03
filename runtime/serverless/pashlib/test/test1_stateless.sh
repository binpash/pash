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

INPUT_FILE="${TMP_DIR}/input_stateless.txt"
OUTPUT_FILE="${TMP_DIR}/output_stateless.txt"
RDV_KEY="ft-test-stateless-$(date +%s)-$RANDOM"

cat > "$INPUT_FILE" <<'EOF'
stateless payload line 1
stateless payload line 2
EOF

"${SCRIPT_DIR}/reader.sh" "$RDV_KEY" "recv-node" "send-node" > "$OUTPUT_FILE" &
READER_PID=$!
sleep 1

cat "$INPUT_FILE" | "${SCRIPT_DIR}/sender.sh" "$RDV_KEY" "send-node" "recv-node" "true" 0 "job-stateless" "folder-stateless" "script-stateless"

wait "$READER_PID"
READER_PID=""

if cmp -s "$INPUT_FILE" "$OUTPUT_FILE"; then
  echo "test1_stateless: PASS"
else
  echo "test1_stateless: FAIL (output mismatch)" >&2
  exit 1
fi
