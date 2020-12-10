#!/bin/bash

PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
WIDTH=0
  
##
# Friendlier wrapper around the pash core (to decouple the two until the project
# stabilizes).
##
usage() {
  cat <<EOF
  Parallelize shell scripts:

    ${0} [-wh] script.sh

    * -w: parallelism width
    * -h: shows this message
EOF
}

# Check which argument we have
while getopts ":w:h" opt; do
  case $opt in
    w)
      WIDTH=${OPTARG}
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
  if [ $WIDTH = 0 ]; then
    echo python3.8 $PASH_TOP/compiler/pash.py $file
  else
    echo python3.8 $PASH_TOP/compiler/pash.py --split_fan_out 2 $file
  fi
done
