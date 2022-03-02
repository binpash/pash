#! /usr/bin/env bash
# Infer Unix-like distribution
#
# TODO: This logic appears verbatim in a few places.
# Deduplicate by having other scripts call this one.

if type lsb_release >/dev/null 2>&1 ; then
    distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
    distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]'
