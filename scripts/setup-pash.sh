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
# remove the folder in case it exists
rm -rf $PYTHON_PKG_DIR
# create the new folder
mkdir -p $PYTHON_PKG_DIR

echo "Installing python dependencies..."

python3 -m pip install -r "$PASH_TOP/requirements.txt" --no-cache-dir --root $PYTHON_PKG_DIR --ignore-installed

## numpy and matplotlib are only needed to generate the evaluation plots so they should not be in the main path
if [[ "$install_eval" == 1 ]];  then
    python3 -m pip install numpy --root $PYTHON_PKG_DIR --ignore-installed #&> $LOG_DIR/pip_install_numpy.log
    python3 -m pip install matplotlib --root $PYTHON_PKG_DIR --ignore-installed #&> $LOG_DIR/pip_install_matplotlib.log
fi

# clean the python packages
cd $PYTHON_PKG_DIR
# can we find a better alternative to that                                      
pkg_path=$(find . \( -name "site-packages" -or -name "dist-packages" \) -type d)
for directory in $pkg_path; do
  # using which to avoid the `-i` alias in many distros
  $(which cp) -r $directory/* ${PYTHON_PKG_DIR}/
done


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

