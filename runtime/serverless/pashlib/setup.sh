#!/bin/bash
set -euo pipefail

curl https://sh.rustup.rs -sSf | sh
source "$HOME/.cargo/env"