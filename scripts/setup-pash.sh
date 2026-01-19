#!/usr/bin/env bash

set -e

## TODO: Maybe hide stdout and stderr to logs by default and only if debug flag exists show

cd "$(dirname "$0")"
# set PASH_TOP
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

cd $PASH_TOP
. "$PASH_TOP/scripts/utils.sh"
read_cmd_args $@

LOG_DIR=$PASH_TOP/install_logs
mkdir -p $LOG_DIR
PYTHON_PKG_DIR=$PASH_TOP/python_pkgs

# Check if venv module is available
if ! python3 -m venv --help &> /dev/null; then
    echo "Error: python3-venv module not found. Please install it:" >&2
    echo "  Ubuntu/Debian: sudo apt install python3-venv" >&2
    echo "  Fedora/RHEL: sudo dnf install python3" >&2
    exit 1
fi


# Create virtual environment
echo "Creating virtual environment..."
python3 -m venv $PYTHON_PKG_DIR

# Install dependencies
echo "Installing Python dependencies..."
$PYTHON_PKG_DIR/bin/pip install --upgrade pip
$PYTHON_PKG_DIR/bin/pip install -r "$PASH_TOP/requirements.txt"

## numpy and matplotlib are only needed to generate the evaluation plots so they should not be in the main path
if [[ "$install_eval" == 1 ]]; then
    $PYTHON_PKG_DIR/bin/pip install numpy matplotlib
fi

# Build runtime tools: eager, split
echo "Building runtime tools..."
cd "$PASH_TOP/runtime/"
case "$distro" in
    freebsd*) 
        gmake #&> $LOG_DIR/make.log
        ;;
    *)
        make #&> $LOG_DIR/make.log
        if [ -f /.dockerenv ]; then
            # issue with docker only
            python3 -m pip install -U --force-reinstall pip
            cp "$PASH_TOP"/pa.sh /usr/bin/
        fi
        ;;
esac

echo "Generating input files..."
$PASH_TOP/evaluation/tests/input/setup.sh

# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
echo " * * * "
echo "Do not forget to export PASH_TOP before using pash: \`export PASH_TOP=$PASH_TOP\`"
echo " * * * "
# in case we are running on docker or CI, installation is complete at this moment
if [[ -f /.dockerenv || -f /.githubenv ]]; then
    exit 0  
fi

