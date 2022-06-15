#! /usr/bin/env bash

# Build, install, and use PaSh in a new Docker container.
#
#   deploy.sh IMAGE VERSION FORMAT
#
# Exit code matches that from running `./tests/$3/main.sh` in the
# container.  0 is a sign of PaSh's functionality in the other
# environment.
#
# Scripts in ./tests may use /package in the container to store
# information on the host, but printing should suffice.
#
# $1 - IMAGE name from Docker Hub to use
# $2 - VERSION of PaSh to install in a container
# $3 - FPM package FORMAT for the named PaSh version

set -euo pipefail

me="$(readlink -f "${BASH_SOURCE[0]}")"
mydirname="$(dirname "$me")"
myname="$(basename "$me")"

if [ -f /.dockerenv ]; then
    set +x
    # Pass control to vertification script.
    target_version="$1"
    output_format="$2"

    # Keep this here, so tests don't need to keep setting the variable.
    export PASH_TOP=/usr/lib/pash

    "$mydirname/fire-hook.sh" "before-test" "$output_format" "$target_version"
 
    # Basic example must always work, no matter where we're going.
    [ -d "$PASH_TOP" ]
    test_script="${PASH_TOP}/evaluation/intro/hello-world.sh"
    set -x
    time pa.sh "$test_script"
    time bash "$test_script"
    set +x
else
    image="$1"
    target_version="$2"
    output_format="$3"

    # Perform a one-off build for the package.
    "$mydirname/shell.sh" "$target_version" "$output_format"

    # It's important to use Docker hub here because we want to verify
    # PaSh's behavior against other people's work, namely official OS
    # images.
    docker pull "$image"

    # Invoke self in the container.
    docker run \
	   --rm \
	   --volume "${mydirname}:/package" \
	   "$image" \
	   "/package/${myname}" "$target_version" "$output_format"
fi
