#!/usr/bin/env bash
set -euo pipefail

IMAGE_BUILDER_DIR="/src/image-builder"
RUNTIME_SOURCE_DIR="/src/runtime"
OUTPUT_DIR="/opt/pash-built/runtime"

mkdir -p "$OUTPUT_DIR"

find "$RUNTIME_SOURCE_DIR" -maxdepth 1 -type f \
  \( -name '*.c' -o -name '*.h' -o -name '*.sh' -o -name 'Makefile' \) \
  -exec cp -Lf {} "$OUTPUT_DIR"/ \;

make -C "$OUTPUT_DIR"
