#
# LICENSE AT END OF FILE
#
# This is a script to automatically install pacaur from the AUR. It is
# intended for fresh systems with no other means to install from AUR.
# Theoretically, this script can install other packages too. Just modify
# the $aurpkgs variable below.
#
#
# TODO makedeb-makepkg needs to be installed
# which packages to install from AUR, in this order!
# aurpkgs="$(cat packages | tr '\n' ' ' )"

# exit on errors
set -e

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
    mkcd "folder/$pkg"

    info "fetch PKGBUILD for $pkg"
    curl --insecure -o  PKGBUILD "$pkgbuild=$pkg" || echo ' '

    #info "fetch required pgp keys from PKGBUILD"
    #gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    info "make and install ..."
    makedeb-makepkg --format-makedeb -d || echo 'failed'
    cd -
}

export -f package_build_aux
# loop over required packages
for pkg in $(cat packages | tr '\n' ' ' ); do
    package_build_aux $pkg > "logs"/"$pkg.log"
done

echo "done"