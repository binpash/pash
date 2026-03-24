#!/bin/bash

sudo apt-get update

sudo apt-get install -y --no-install-recommends \
  tcpdump curl wget coreutils diffutils gzip bcftools gawk unzip git \
  jq \
  coreutils \
  gawk \
  cmake \
  build-essential \
  libjansson-dev \
  libpcap-dev \
  tar \
  git \
  python3 \
  q-text-as-data \
  grep \
  sed

# Set GO_VERSION *before* using it
GO_VERSION=1.24.2
echo "Installing Go $GO_VERSION"

TOP=$(git rev-parse --show-toplevel)
eval_dir="$TOP/evaluation/benchmarks/analytics"
go_install_dir="${eval_dir}/go_install"

mkdir -p "$go_install_dir"
curl -LO "https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz"
tar -C "$go_install_dir" -xzf "go${GO_VERSION}.linux-amd64.tar.gz"
rm -f "go${GO_VERSION}.linux-amd64.tar.gz"

export GOROOT="$go_install_dir/go"
export GOPATH="$HOME/go"
export PATH="$GOROOT/bin:$GOPATH/bin:$PATH"

# Confirm Go is now working
go version || { echo "Go installation failed"; exit 1; }

# Install zannotate
go install github.com/zmap/zannotate/cmd/zannotate@latest

# Confirm zannotate is now on PATH
command -v zannotate || { echo "zannotate not found on PATH after go install"; exit 1; }
