#! /usr/bin/env bash
# Source this in any new ~/pash-init.sh

set -u
export distro="$("$PASH_TOP/scripts/distro.sh")"

# Adapt to Docker
if [ -f /.dockerenv ]; then
    export LC_ALL='en_US.UTF-8'
    export LANG='en_US.UTF-8'
    export LANGUAGE='en_US.UTF-8'
fi
