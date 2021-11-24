#!/bin/bash
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input/packages}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/input/output/packages}
LOGS=${OUT}/logs
mkdir -p ${OUT} ${LOGS}

info() { echo -e "\e[1m--> $@\e[0m"; }
mkcd() { mkdir -p "$1" && cd "$1"; }

mkdir -p logs
# check if not running as root
# test "$UID" -gt 0 || { info "don't run this as root!"; exit; }

# set link to plaintext PKGBUILDs
pkgbuild="https://aur.archlinux.org/cgit/aur.git/plain/PKGBUILD?h"
info "using '$pkgbuild=<package>' for plaßintext PKGBUILDs"

package_build_aux() {
    pgk=$1
    info "create subdirectory for $pkg"
    mkcd "${OUT}/$pkg"

    info "fetch PKGBUILD for $pkg"
    curl --insecure -o  PKGBUILD "$pkgbuild=$pkg" || echo ' '

    #info "fetch required pgp keys from PKGBUILD"
    #gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    info "make and install ..."
    #makedeb-makepkg --format-makedeb -d || echo 'failed'
    cd -
}

export -f package_build_aux
# loop over required packages
for pkg in $(cat ${IN} | tr '\n' ' ' ); do
    echo $pkg
    package_build_aux $pkg > "${LOGS}"/"$pkg.log"
done
