#!/usr/bin/env bash

distro=${1??Distro not given}
# now do different things depending on distro
case "$distro" in
    freebsd*)  
        # bsd is so fun :)
        tmp=$(TMPDIR=$PASH_TMP_PREFIX mktemp -t pash_XXXXXXXXXX)
        echo "${tmp}"
        ;;
    *)
        mktemp --tmpdir="$PASH_TMP_PREFIX" -u pash_XXXXXXXXXX
        ;;
esac
