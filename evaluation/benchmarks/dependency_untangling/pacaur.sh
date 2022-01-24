#!/bin/bash
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/packages}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/output/packages}
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
    curl --insecure -o  PKGBUILD "$pkgbuild=$pkg" 2> /dev/null|| echo ' '

    #info "fetch required pgp keys from PKGBUILD"
    #gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    info "make and install ..."
    timeout 100 makedeb-makepkg --format-makedeb -d 2>/dev/null|| echo 'failed'
    cd -
}

export -f run_tests
pkg_count=0
# loop over required packages
for pkg in $(cat ${IN} | tr '\n' ' ' ); 
do
    pkg_count=$((pkg_count + 1))
    run_tests $pkg > "${LOGS}"/"$pkg_count.log"
done

echo 'done';
