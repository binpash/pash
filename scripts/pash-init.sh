#! /usr/bin/env bash
# Source this in any new ~/pash-init.sh

set -u
. "$PASH_TOP/scripts/utils.sh"

export distro="$(infer_unix_like_distro)"

# Adapt to Docker
if [ -f /.dockerenv ]; then
    export LC_ALL='en_US.UTF-8'
    export LANG='en_US.UTF-8'
    export LANGUAGE='en_US.UTF-8'
fi
