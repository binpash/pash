#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LEASH_ENTRIES=1000 exec "${SCRIPT_DIR}/../run-leash-benchmark-matrix.sh" --benchmark nlp --parallel_pipelines --parallel_pipelines_limit 1024 "$@"