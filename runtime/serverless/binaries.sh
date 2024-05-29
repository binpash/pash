#!/bin/bash

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

cd "$(dirname "$0")" || exit

BINS=(
    dgsh-tee
    r_merge
    r_split
    r_unwrap
    r_wrap
    wait_for_output_and_sigpipe_rest.sh
)

SRC_DIR=$PASH_TOP/runtime
DST_DIR=$PASH_TOP/runtime/serverless/runtime

for BIN in "${BINS[@]}"; do
    cp "$SRC_DIR"/"$BIN" "$DST_DIR"/"$BIN"
done
