#!/usr/bin/env bash
if type lsb_release >/dev/null 2>&1 ; then
    distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
    distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

# convert to lowercase
distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# now do different things depending on distro
case "$distro" in
    freebsd*)  
        # bsd is so fun :)
        mkdir -p "$PASH_TMP_PREFIX"
        tmp=$(mktemp -u pash_XXXXXXXXXX)
        mv $tmp "$PASH_TMP_PREFIX"
        echo "${PASH_TMP_PREFIX}/${tmp}"
        ;;
    *)
        echo "$(mktemp --tmpdir="$PASH_TMP_PREFIX" -u pash_XXXXXXXXXX)"
        ;;
esac
