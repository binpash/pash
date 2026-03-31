#!/usr/bin/env bash
set -euo pipefail

RUNTIME_DIR="/var/task/runtime"
RUNTIME_LIB_DIR="$RUNTIME_DIR/lib"

export PATH="$PATH:$RUNTIME_DIR"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:$RUNTIME_LIB_DIR"

check_binary() {
  local name=$1
  shift

  echo "[test.sh] checking $name"
  command -v "$name" >/dev/null
  "$@"
}

check_present() {
  local name=$1
  echo "[test.sh] checking $name"
  local path
  path="$(command -v "$name")"
  [[ -n "$path" ]]
  [[ -x "$path" ]]
}

check_runtime_path() {
  local relative_path=$1
  echo "[test.sh] checking $relative_path"
  [[ -x "/var/task/$relative_path" ]]
}

check_absent() {
  local name=$1
  echo "[test.sh] checking $name is absent"
  ! command -v "$name" >/dev/null 2>&1
}

check_binary openssl openssl version
check_binary ffmpeg ffmpeg -version
check_binary gzip gzip --version
check_binary convert convert --version
check_binary file file --version
check_binary xargs xargs --version
check_binary tcpdump tcpdump --version
check_binary col col --help
check_binary rev rev --help
check_binary which which --version

check_runtime_path runtime/eager
check_runtime_path runtime/split
check_runtime_path runtime/r_merge
check_runtime_path runtime/r_split
check_runtime_path runtime/r_wrap
check_runtime_path runtime/r_unwrap
check_runtime_path runtime/set-diff
check_runtime_path runtime/dgsh-tee
check_runtime_path runtime/pashlib
check_runtime_path runtime/pashlib-ft

echo "[test.sh] exercising convert JPEG coder"
convert -size 1x1 xc:white /tmp/test-image.jpg
file /tmp/test-image.jpg | grep -qi 'jpeg'

# The builder stage may install toolchain packages, but the final runtime image
# should not carry them.
check_absent gcc
check_absent git
check_absent make
