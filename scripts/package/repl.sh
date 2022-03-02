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

cd "$(dirname "${BASH_SOURCE[0]}")" # Help limit filesystem writes to this directory.
mkdir -p output

main() {
    docker run \
	   --entrypoint /bin/bash \
	   --interactive \
	   --tty \
	   --rm \
	   --user "$(id -u):$(id -g)" \
	   --volume "$PASH_TOP/scripts/package/output:/out" \
	   --volume "$PASH_TOP/scripts/package/tools:/tools" \
	   --volume "$PASH_TOP:/src" \
	   fpm --rcfile /tools/start -i
}

# Only enter if script was not sourced
(return 2>/dev/null) || main "$@"
