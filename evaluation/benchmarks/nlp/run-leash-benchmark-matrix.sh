#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/../run-leash-benchmark-matrix.sh" --benchmark nlp "$@"
