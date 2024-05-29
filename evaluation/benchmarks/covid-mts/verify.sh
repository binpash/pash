#!/bin/bash

map_outfile_to_target_hashfile() {
    # Remove the first "__", then remove until the last "__",
    # then remove until the first "."
    # then change ".out" to ".hash"
    echo "$1" |
        sed "s/__.*__[^.]*//" |
        sed "s/\.out/\.hash/"
}

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

[ -z "$AWS_ACCOUNT_ID" ] && {
  echo "AWS_ACCOUNT_ID not set"
  exit
}

[ -z "$AWS_BUCKET" ] && {
  echo "AWS_BUCKET not set"
  exit
}

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/covid-mts"
OUTPUTS_DIR="$BENCHMARK_DIR/outputs"
HASHES_DIR="$BENCHMARK_DIR/hashes"

mkdir -p hashes/

if [[ "$*" == *"--generate"* ]]
then
    for GROUND_TRUTH_PATH in "$OUTPUTS_DIR"/*Local*Bash*.out
    do
        GROUND_TRUTH=$(basename "$GROUND_TRUTH_PATH")
        HASH_FILE=$(map_outfile_to_target_hashfile "$GROUND_TRUTH")
        HASH_PATH="$HASHES_DIR/$HASH_FILE"

        HASH=$(sha256sum "$GROUND_TRUTH_PATH" | cut -d" " -f1)

        echo Writing to "$HASH_PATH"

        echo "$HASH" >"$HASH_PATH"
    done
else
    for OUTPUT_PATH in "$OUTPUTS_DIR"/*.out
    do
        OUTPUT=$(basename "$OUTPUT_PATH")
        TARGET_HASH_FILE=$(map_outfile_to_target_hashfile "$OUTPUT")
        TARGET_HASH_PATH="$HASHES_DIR/$TARGET_HASH_FILE"
        HASH_FILE=${OUTPUT%.*}.hash
        HASH_PATH="$HASHES_DIR/$HASH_FILE"

        HASH=$(sha256sum "$OUTPUT_PATH" | cut -d" " -f1)

        echo "$HASH" >"$HASH_PATH"

        if diff "$HASH_PATH" "$TARGET_HASH_PATH" >/dev/null
        then
            echo "Hashes for $OUTPUT match"
        else
            echo "Hashes for $OUTPUT do not match"
        fi
    done
fi
