#!/usr/bin/env bash

pkgs="p7zip-full curl wget unzip"

sudo apt-get update
for pkg in $pkgs; do
    if ! dpkg -s "$pkg" > /dev/null 2>&1; then
        sudo apt-get install -y --no-install-recommends "$pkg"
    fi
done

# Install pandoc if not installed
if ! dpkg -s pandoc > /dev/null 2>&1; then
    # pandoc v2.x does not support arm64; use v3.5
    arch=$(dpkg --print-architecture)
    wget https://github.com/jgm/pandoc/releases/download/3.5/pandoc-3.5-1-"${arch}".deb
    sudo dpkg -i pandoc-3.5-1-"${arch}".deb || sudo apt-get install -f -y --no-install-recommends
    rm pandoc-3.5-1-"${arch}".deb
fi

# Install Node.js (18.x) and npm via NodeSource if not present
if ! command -v node > /dev/null 2>&1; then
    curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -
    sudo apt-get install -y --no-install-recommends nodejs
fi

if ! command -v node > /dev/null 2>&1; then
    echo "Node.js installation failed." >&2
    exit 1
fi

# Install npm dependencies for stem-words.js
input_dir="$(dirname "$0")/input"
cd "$input_dir" || exit 1
if [ ! -d node_modules ]; then
    npm install
fi

# Download inputs (small dataset)
benchmark_dir="$(dirname "$0")"
cd "$benchmark_dir" || exit 1
./inputs.sh --small
