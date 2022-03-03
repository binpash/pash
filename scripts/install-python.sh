#! /usr/bin/env bash

set -euo pipefail
PYTHON_PKG_DIR="$PASH_TOP/python_pkg_root"
python3 -m pip install --root "$PYTHON_PKG_DIR" --ignore-installed -r "$PASH_TOP/requirements.txt"
ln -frs "$(find "$PYTHON_PKG_DIR" -type d -name site-packages)" "$PASH_TOP/python_pkgs"
