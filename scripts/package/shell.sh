#! /usr/bin/env bash
# Launch REPL for building PaSh system packages.

here="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
PASH_TOP=$(readlink -f "$here/../..")

# We need FPM's source code to build a key Docker image.
printf "COMPARE\n  >$here\n  >%s\n " "$(pwd)"
if [ ! -d "$here/fpm" ]; then
    git clone --depth 1 git@github.com:jordansissel/fpm.git "$here/fpm"
fi

# Make a Docker image that knows about building packages.
if ! docker image inspect fpm 2>&1 >/dev/null; then
    make -C fpm docker-release-everything
fi

mkdir -p "$here/output"

main() {
    local docker_interactive_flags=''
    local entrypoint_interactive_flags=''

    if [[ $- == *i* ]]; then
        set -x
        docker_interactive_flags='--interactive --tty'
        entrypoint_interactive_flags='--rcfile'
        set +x
    fi

    docker run \
	   --entrypoint /bin/bash \
	   --rm \
	   --user "$(id -u):$(id -g)" \
	   --volume "$here/output:/out" \
	   --volume "$here/tools:/tools" \
	   --volume "$PASH_TOP:/src" \
           $docker_interactive_flags \
           fpm \
           $entrypoint_interactive_flags \
           /tools/start "$@"
}

# Only enter if script was not sourced
(return 2>/dev/null) || main "$@"
