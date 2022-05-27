#!/usr/bin/env bash

distro=${1??Distro not given}
# echo "$PASH_TMP_PREFIX/pash_$RANDOM$RANDOM$RANDOM"
mktemp -u "$PASH_TMP_PREFIX/pash_XXXXXXXXXX"
