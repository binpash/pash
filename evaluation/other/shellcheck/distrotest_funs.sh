#!/bin/bash

die()
{
    echo "$*" >&4;
    exit 1;
}

distrotest_loop()
{
    while read -r distro setup
    do
        [[ "$distro" = "#"* || -z "$distro" ]] && continue

        printf '%s ' "$distro" # >&3
        docker pull "$distro" || die "Can't pull $distro"
        printf 'pulled. ' # >&3

        tmp=$(mktemp -d) || die "Can't make temp dir"
        cp -r "${SHELLCHECK_DIR}" "$tmp/" || die "Can't populate test dir"
        printf 'Result: ' # >&3
        < /dev/null docker run -v "$tmp:/mnt" "$distro" sh -c "
        $setup
        cd /mnt/shellcheck || exit 1
        test/buildtest
        "
        ret=$?
        if [ "$ret" = 0 ]
        then
            echo "OK" # >&3
        else
            echo "FAIL with $ret. See $log" # >&3
            final=1
        fi
        rm -rf "$tmp"
    done
}

export -f die
export -f distrotest_loop
