#!/usr/bin/env bash

set -e

cd "$(dirname "$0")"
SCRIPT_DIR="$(pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

. "$SCRIPT_DIR/utils.sh"
read_cmd_args "$@"

LOG_DIR=$REPO_ROOT/install_logs
mkdir -p "$LOG_DIR"
PYTHON_PKG_DIR=$REPO_ROOT/python_pkgs

# Check if venv module is available
if ! python3 -m venv --help &> /dev/null; then
    echo "Error: python3-venv module not found. Please install it:" >&2
    echo "  Ubuntu/Debian: sudo apt install python3-venv" >&2
    echo "  Fedora/RHEL: sudo dnf install python3" >&2
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d "$PYTHON_PKG_DIR" ]; then
    echo "Creating virtual environment..."
    python3 -m venv "$PYTHON_PKG_DIR"
fi

# Upgrade pip
echo "Upgrading pip..."
"$PYTHON_PKG_DIR/bin/pip" install --upgrade pip

# Install PaSh in editable mode (this also builds runtime binaries via setup.py)
echo "Installing PaSh and dependencies..."
"$PYTHON_PKG_DIR/bin/pip" install -e "$REPO_ROOT"

# Install evaluation dependencies if requested
if [[ "$install_eval" == 1 ]]; then
    echo "Installing evaluation dependencies..."
    "$PYTHON_PKG_DIR/bin/pip" install numpy matplotlib
fi

# Generate input files for tests
echo "Generating input files..."
"$REPO_ROOT/evaluation/tests/input/setup.sh"

# Docker-specific setup
if [ -f /.dockerenv ]; then
    # Make pash available system-wide in Docker
    cp "$REPO_ROOT/pa.sh" /usr/bin/ 2>/dev/null || true
fi

echo " * * * "
echo "PaSh installation complete!"
echo ""
echo "To use PaSh, either:"
echo "  1. Activate the virtual environment: source $PYTHON_PKG_DIR/bin/activate"
echo "     Then run: pash <script.sh>"
echo ""
echo "  2. Or use pa.sh directly: $REPO_ROOT/pa.sh <script.sh>"
echo ""
echo "  3. Or set PASH_TOP and add to PATH:"
echo "     export PASH_TOP=$REPO_ROOT/src/pash"
echo "     export PATH=\$PATH:$REPO_ROOT"
echo " * * * "

# In CI/Docker environments, exit without prompting
if [[ -f /.dockerenv || -f /.githubenv ]]; then
    exit 0
fi
