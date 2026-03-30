#!/usr/bin/env bash
set -euo pipefail

IMAGE_BUILDER_DIR="${IMAGE_BUILDER_DIR:-/opt/image-builder}"

read_list_file() {
  local file_path=$1
  [[ -f "$file_path" ]] || return 0

  awk '
    {
      sub(/[[:space:]]*#.*/, "", $0)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
      if (length($0) > 0) {
        print $0
      }
    }
  ' "$file_path"
}

binary_name_from_entry() {
  local entry=$1
  printf '%s\n' "${entry%%=*}"
}

LOCAL_RUNTIME_BINARIES=(
  openssl
  ffmpeg
  gzip
  convert
  file
  xargs
  tcpdump
  col
  rev
  which
  eager
  split
  r_merge
  r_split
  r_wrap
  r_unwrap
  set-diff
  dgsh-tee
  pashlib
  pashlib-ft
)

declare -A SEEN=()
BINARIES_TO_CHECK=()

for binary_name in "${LOCAL_RUNTIME_BINARIES[@]}"; do
  if [[ -z "${SEEN[$binary_name]:-}" ]]; then
    BINARIES_TO_CHECK+=("$binary_name")
    SEEN["$binary_name"]=1
  fi
done

while IFS= read -r entry; do
  binary_name="$(binary_name_from_entry "$entry")"
  if [[ -z "${SEEN[$binary_name]:-}" ]]; then
    BINARIES_TO_CHECK+=("$binary_name")
    SEEN["$binary_name"]=1
  fi
done < <(read_list_file "$IMAGE_BUILDER_DIR/binaries.txt")

for binary_name in "${BINARIES_TO_CHECK[@]}"; do
  runtime_binary_path="runtime/$binary_name"
  if [[ ! -x "$runtime_binary_path" ]]; then
    echo "Missing runtime binary: $runtime_binary_path" >&2
    exit 1
  fi

  echo "[verify-runtime-linking] checking $runtime_binary_path"
  ldd_output="$(ldd "$runtime_binary_path" 2>&1 || true)"

  if grep -q "not found" <<<"$ldd_output"; then
    echo "$ldd_output" >&2
    echo "Unresolved shared library dependency for $runtime_binary_path" >&2
    exit 1
  fi
done

echo "[verify-runtime-linking] exercising convert JPEG coder"
convert -size 1x1 xc:white /tmp/verify-runtime-linking.jpg
