#! /usr/bin/env bash
# Launch REPL for building PaSh system packages.

# We need FPM's source code to build a key Docker image.
if [ ! -d fpm ]; then
    git clone --depth 1 git@github.com:jordansissel/fpm.git
fi

# Make a Docker image that knows about building packages.
if ! docker image inspect fpm 2>&1 >/dev/null; then
    make -C fpm docker-release-everything
fi

here="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
PASH_TOP=$(readlink -f "$here/../..")
mkdir -p "$here/output"
echo $here

main() {
    docker run \
	   --entrypoint /bin/bash \
	   --interactive \
	   --tty \
	   --rm \
	   --user "$(id -u):$(id -g)" \
	   --volume "$here/output:/out" \
	   --volume "$here/tools:/tools" \
	   --volume "$PASH_TOP:/src" \
	   fpm --rcfile /tools/start "$@"
}

# Only enter if script was not sourced
(return 2>/dev/null) || main "$@"
