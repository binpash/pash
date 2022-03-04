#! /usr/bin/env bash
# Called by system package manager to request on-site setup on host.

set -euo pipefail

pashd=/usr/lib/pash
logd=/var/log/pash
export PASH_TOP="$pashd"

say_task() {
    printf 'pash: %s (log: %s/%s)\n' "$1" "$logd" "$2"
}

mkdir -vp "$logd"

say_task "python dependencies" "install-python.log"

if ! command -v pip 2>&1 >/dev/null; then
    /bin/python3 -m ensurepip --upgrade
fi

"$pashd/scripts/install-python.sh" 2>&1 | tee "$logd/install-python.log"

cd "$pashd/compiler/parser"
ln -fs "$pashd/pa.sh" /usr/bin/pa.sh

say_task "libdash" "libdash.log"
rm -rf libdash
git clone https://github.com/angelhof/libdash/
make libdash 2>&1 | tee "$logd/libdash.log"

say_task "runtime" "runtime.log"
cd "$pashd/runtime"
make clean
make 2>&1 | tee "$logd/runtime.log"
