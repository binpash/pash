#!/usr/bin/env bash
set -euo pipefail

IMAGE_BUILDER_DIR="/src/image-builder"
RUNTIME_SOURCE_DIR="/src/runtime"
SERVERLESS_RUNTIME_DIR="/src/serverless-runtime"
PASHLIB_SOURCE_DIR="/src/serverless-pashlib"
OUTPUT_DIR="/opt/pash-built/runtime"

mkdir -p "$OUTPUT_DIR"

find "$RUNTIME_SOURCE_DIR" -maxdepth 1 -type f \
  \( -name '*.c' -o -name '*.h' -o -name '*.sh' -o -name 'Makefile' \) \
  -exec cp -Lf {} "$OUTPUT_DIR"/ \;

make -C "$OUTPUT_DIR"

cp -a "$SERVERLESS_RUNTIME_DIR/." "$OUTPUT_DIR/"

# Prefer package-managed userland tools from the final image. The checked-in
# serverless runtime bundle still contains older copies plus support libraries,
# and those can break the RPM-installed binaries once LD_LIBRARY_PATH points at
# /var/task/runtime/lib.
rm -f \
  "$OUTPUT_DIR/col" \
  "$OUTPUT_DIR/convert" \
  "$OUTPUT_DIR/file" \
  "$OUTPUT_DIR/gzip" \
  "$OUTPUT_DIR/openssl" \
  "$OUTPUT_DIR/rev" \
  "$OUTPUT_DIR/tcpdump" \
  "$OUTPUT_DIR/which" \
  "$OUTPUT_DIR/xargs"
rm -rf "$OUTPUT_DIR/lib"

mkdir -p /src/serverless/runtime
cp -a "$OUTPUT_DIR/." /src/serverless/runtime/

(
  cd "$PASHLIB_SOURCE_DIR"
  cargo build --release --bin pashlib --bin pashlib-ft
)

cp -Lf "$PASHLIB_SOURCE_DIR/target/release/pashlib" "$OUTPUT_DIR/pashlib"
cp -Lf "$PASHLIB_SOURCE_DIR/target/release/pashlib-ft" "$OUTPUT_DIR/pashlib-ft"
