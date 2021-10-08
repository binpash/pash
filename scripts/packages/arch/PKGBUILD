# Maintainer: Dimitris Karnikis <dkarnikis@gmail.com>
# Contributor: Fabian Bornschein <fabiscafe-cat-mailbox-dog-org>

pkgbase="pash-shell"
pkgname=('pash-shell' 'pash-shell-docs')
_reponame="pash"
pkgver=4
pkgrel=2
pkgdesc="PaSh: Light-touch Data-Parallel Shell Processing"
url="https://binpa.sh/"
arch=("x86_64")
license=('MIT')
_commit=2b4db9ddb1616ce6fe56575724c966177bb950f7  # tags/v.04^0
source=("git+https://github.com/binpash/pash.git#commit=$_commit")
depends=('python-jsonpickle' 'python-yaml' 'python-numpy' 'python-matplotlib')
makedepends=('git' 'gcc' 'gcc10')
checkdepends=('inetutils')
sha256sums=('SKIP')

pkgver() {
  cd "${srcdir}/${_reponame}"
  git describe --tags | sed 's/v.//;s/^0*//;s/-/+/g'
}

prepare() {
  # Hardcode PASH_TOP, as we already know where it is placed in the FS
  # and there is no need for the user the set it later on
  sed -i 's|^export PASH_TOP=.*|export PASH_TOP="/opt/pash"|' "${srcdir}/${_reponame}/pa.sh"
}

build() {
  export PASH_TOP="${srcdir}/${_reponame}"

  cd ${PASH_TOP}
  git submodule init
  git submodule update

  cd compiler/parser
  make libdash

  cd "${PASH_TOP}/runtime"
  make
}

check() {
  cd "${PASH_TOP}"
  "${PASH_TOP}/evaluation/tests/input/setup.sh"
}

package_pash-shell() {
  cd "${srcdir}/${_reponame}"
  install -Dm644 LICENSE "${pkgdir}/usr/share/licenses/pash-shell/LICENSE"
  /usr/bin/install -dm755 "${pkgdir}/opt/pash"

  mv ./annotations ./compiler ./evaluation ./runtime "${pkgdir}/opt/pash/"
  /usr/bin/install -m755 ./pa.sh "${pkgdir}/opt/pash/pa.sh"
  
  /usr/bin/install -dm755 "${pkgdir}/usr/bin"
  ln -s /opt/pash/pa.sh "${pkgdir}/usr/bin/pa.sh"
}

package_pash-shell-docs() {
  pkgdesc="Documentation for PaSh"
  arch=("any")
  depends=()
  cd "${srcdir}/${_reponame}"
  /usr/bin/install -dm755 "${pkgdir}/opt/pash"
  mv ./docs "${pkgdir}/opt/pash/"
  install -Dm644 LICENSE "${pkgdir}/usr/share/licenses/pash-shell-docs/LICENSE"
}
