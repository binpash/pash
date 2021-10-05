#!/usr/bin/env bash

distro=${1??Distro not given}
TMPDIR=$PASH_TMP_PREFIX mktemp -u -t pash_XXXXXXXXXX
