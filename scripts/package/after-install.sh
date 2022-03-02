#! /usr/bin/env bash
# Called by system package manager to request on-site setup on host.

set -euo pipefail

pashd=/usr/lib/pash
logd=/var/log/pash

say_task() {
    printf 'pash: %s (log: %s/%s)\n' "$1" "$logd" "$2"
}

mkdir -vp "$logd"
cd "$pashd/compiler/parser"
ln -fs "$pashd/pa.sh" /usr/bin/pa.sh
