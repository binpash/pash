#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# manual fix to bypass the compilation err
cd $PASH_TOP
rm {1,2,3}grams 2> /dev/null || true
mkfifo {1,2,3}grams

exec "${SCRIPT_DIR}/../run-leash-benchmark-matrix.sh" --benchmark web-search "$@"

cd $PASH_TOP
rm {1,2,3}grams 2> /dev/null || true