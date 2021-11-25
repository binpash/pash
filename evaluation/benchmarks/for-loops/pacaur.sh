#!/bin/bash
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input/packages}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/input/output/packages}
LOGS=${OUT}/logs
mkdir -p ${OUT} ${LOGS}

info() { echo -e "\e[1m--> $@\e[0m"; }
mkcd() { mkdir -p "$1" && cd "$1"; }

# check if not running as root
# test "$UID" -gt 0 || { info "don't run this as root!"; exit; }

# set link to plaintext PKGBUILDs
pkgbuild="https://aur.archlinux.org/cgit/aur.git/plain/PKGBUILD?h"

run_tests() {
    pgk=$1
    info "create subdirectory for $pkg"
    mkcd "${OUT}/$pkg"

    info "fetch PKGBUILD for $pkg"
    curl --insecure -o  PKGBUILD "$pkgbuild=$pkg" || echo ' '

    #info "fetch required pgp keys from PKGBUILD"
    #gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    info "make and install ..."
    timeout 100 makedeb-makepkg --format-makedeb -d || echo 'failed'
    cd -
}

export -f run_tests
# loop over required packages
for pkg in $(cat ${IN} | tr '\n' ' ' ); do
    run_tests $pkg > "${LOGS}"/"$pkg.log" 2> /dev/null
done

echo 'done';
