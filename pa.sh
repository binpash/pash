#!/bin/bash

export PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
export PASH_PARSER=${PASH_TOP}/parser/parse_to_json.native
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

WIDTH=2
VIEW_ONLY=""
  
##
# Friendlier wrapper around the pash core (to decouple the two until the project
# stabilizes).
##
usage() {
  cat <<EOF
Parallelize shell scripts:

	  ${0} [-whv] script.sh

* -h, --help         show this help message
* -w n, --width n    configure the level of parallelism sought
* -v, --view-only    only view parallel script, not execute
EOF
}

# Check which argument we have
while getopts ":w:hv" opt; do
  case $opt in
    w)
      WIDTH=${OPTARG}
      ;;
    v)
      VIEW_ONLY="--compile_optimize_only"
      ;;
    h)
      usage;
      exit 0;
      ;;
    *)
      usage;
      exit 1;
      ;;
  esac
done

shift $(( OPTIND - 1 ))

if [ $# = 0 ]; then
  echo "Interactive mode not supported yet."
  exit 1;
fi

for file in "$@"; do
    python3 $PASH_TOP/compiler/pash.py --log_file /tmp/pash.log --split_fan_out $WIDTH $VIEW_ONLY $file
done
