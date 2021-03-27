#!/bin/bash
# This script runs 'buildtest' on each of several distros
# via Docker.
set -o pipefail

exec 3>&1 4>&2

[ -e "${SHELLCHECK_DIR}/ShellCheck.cabal" ] || die "ShellCheck.cabal not in this dir"

# [ "$1" = "--run" ] || {
# cat << EOF
# This script pulls multiple distros via Docker and compiles
# ShellCheck and dependencies for each one. It takes hours,
# and is still highly experimental.

# Make sure you're plugged in and have screen/tmux in place,
# then re-run with $0 --run to continue.

# Also note that dist* will be deleted.
# EOF
# exit 0
# }

echo "Deleting 'dist' and 'dist-newstyle'..."
rm -rf dist dist-newstyle

log=$(mktemp) || die "Can't create temp file"
date >> "$log" || die "Can't write to log"

echo "Logging to $log" >&3
## If I keep this on, the script output (together with Dish output) is
## redirected
# exec >> "$log" 2>&1

final=0

cat $IN | distrotest_loop

# distrotest_loop << EOF
# # Docker tag          Setup command
# debian:stable         apt-get update && apt-get install -y cabal-install
# debian:testing        apt-get update && apt-get install -y cabal-install
# ubuntu:latest         apt-get update && apt-get install -y cabal-install
# haskell:latest        true
# opensuse/leap:latest  zypper install -y cabal-install ghc
# fedora:latest         dnf install -y cabal-install ghc-template-haskell-devel findutils
# archlinux/base:latest pacman -S -y --noconfirm cabal-install ghc-static base-devel

# # Other versions we want to support
# ubuntu:18.04          apt-get update && apt-get install -y cabal-install

# # Misc Haskell including current and latest Stack build
# ubuntu:18.04          set -e; apt-get update && apt-get install -y curl && curl -sSL https://get.haskellstack.org/ | sh -s - -f && cd /mnt && exec test/stacktest
# EOF

exit "$final"
