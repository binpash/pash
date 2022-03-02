#! /usr/bin/env bash
# Source this in any new ~/pash-init.sh

export distro="$(/usr/local/pash/scripts/distro.sh)"

# Adapt to Docker
if [ -f /.dockerenv ]; then
    export PASH_TOP=/opt/pash
    export LC_ALL='en_US.UTF-8'
    export LANG='en_US.UTF-8'
    export LANGUAGE='en_US.UTF-8'
fi
