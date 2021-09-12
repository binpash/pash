#!/usr/bin/env bash

export PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

if [ "$#" -eq 1 ] && [ "$1" = "--init" ]; then
  $PASH_TOP/compiler/superoptimize.sh
  exit
fi

if ! command -v python3 &> /dev/null
then
    echo "Python >=3 could not be found"
    exit
fi

mkfifo /tmp/fifo

# TODO: Do we need setsid?
{ tee /tmp/pipo > /tmp/fifo <&3 3<&- & } 3<&0

PASH_FROM_SH="pa.sh" python3 $PASH_TOP/compiler/pash.py "$@" < /tmp/fifo
# TODO: Does interactivity mean more than just reading line by line? Do we need to also return something back?

rm /tmp/fifo