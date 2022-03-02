#! /usr/bin/env bash
# Called by system package manager to request on-site setup on host.

set -euo pipefail

pashd=/usr/lib/pash
logd=/var/log/pash

ln -s "$pashd/pa.sh" /usr/bin/pa.sh

say() {
    printf 'pash: after-install: '
    printf $@
    printf '\n'
}

say_task() {
    say '%s (log: %s/%s)' "$1" "$logd" "$2"
}

mkdir -p "$logd"
cd "$pashd/compiler/parser"

case $(bash "$pashd/scripts/distro.sh") in
    freebsd*)
	cp Makefile Makefile.backup
        gsed -i 's/ make/ gmake/g' Makefile

	say_task "compile libdash (FreeBSD)" libdash.log
        gmake libdash &> "$logd/libdash.log"

	say_task "compile runtime (FreeBSD)" runtime.log
        gmake -C ../../runtime/ &> "$logd/runtime.log"
        ;;
    *)
	say_task "compile libdash" libdash.log
        make libdash &>"$logd/libdash.log"

	say_task "compile runtime" runtime.log
        make -C ../../runtime/ &>"$logd/runtime.log"
        ;;
esac

wait
