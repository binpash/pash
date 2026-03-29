#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
SERVERLESS_DIR="$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)"
PASH_TOP="$(CDPATH= cd -- "$SERVERLESS_DIR/../.." && pwd)"
IMAGE_TAG="${1:-pash-serverless:latest}"
CONTEXT_DIR="$(mktemp -d "$SERVERLESS_DIR/.docker-build-context.XXXXXX")"

cleanup() {
  rm -rf "$CONTEXT_DIR"
}
trap cleanup EXIT

mkdir -p "$CONTEXT_DIR/runtime" "$CONTEXT_DIR/serverless" "$CONTEXT_DIR/image-builder" "$CONTEXT_DIR/vendored-runtime"

while IFS= read -r source_path; do
  cp -Lf "$source_path" "$CONTEXT_DIR/runtime/"
done < <(find "$PASH_TOP/runtime" -maxdepth 1 -type f \
  \( -name '*.c' -o -name '*.h' -o -name '*.sh' -o -name 'Makefile' \) | sort)

cp -a "$SERVERLESS_DIR/aws" "$CONTEXT_DIR/serverless/aws"
cp -Lf "$SERVERLESS_DIR/lambda-function.py" "$CONTEXT_DIR/serverless/lambda-function.py"
cp -Lf "$SERVERLESS_DIR/test.sh" "$CONTEXT_DIR/serverless/test.sh"
cp -Lf "$SERVERLESS_DIR/runtime/ffmpeg" "$CONTEXT_DIR/vendored-runtime/ffmpeg"
cp -a "$SCRIPT_DIR/." "$CONTEXT_DIR/image-builder/"

docker build -f "$CONTEXT_DIR/image-builder/Dockerfile" -t "$IMAGE_TAG" "$CONTEXT_DIR"

# Fail the build early if any packaged binary still has unresolved shared libs.
docker run --rm --entrypoint /bin/bash "$IMAGE_TAG" \
  /opt/image-builder/verify-runtime-linking.sh
