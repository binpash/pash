#!/usr/bin/env bash

distro=${1??Distro not given}
# now do different things depending on distro
tmp_file=$(TMPDIR=$PASH_TMP_PREFIX mktemp -t pash_XXXXXXXXXX)
echo "${tmp_file}"
