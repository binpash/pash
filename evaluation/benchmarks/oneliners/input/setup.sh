#!/bin/bash

#set -e

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

# another solution for capturing HTTP status code
# https://superuser.com/a/590170
input_files="1M.txt 10M.txt 100M.txt 1G.txt dict.txt 3G.txt 10G.txt 100G.txt all_cmds.txt all_cmdsx100.txt small"

if [[ "$1" == "-c" ]]; then
    rm -rf $input_files
    exit
fi

setup_dataset() {
  if [ "$#" -eq 1 ] && [ "$1" = "--small" ]; then
    if [ ! -d ./small ]; then
      echo "Generating small-size inputs"
      # FIXME PR: Do we need all of them?
      curl -sf 'http://pac-n4.csail.mit.edu:81/pash_data/small/oneliners.zip' > oneliners.zip
      unzip oneliners.zip
      rm -f oneliners.zip
    fi
    return 0
  fi

    if [ ! -f ./1M.txt ]; then
        curl -sf 'http://ndr.md/data/dummy/1M.txt' > 1M.txt
        if [ $? -ne 0 ]; then
            echo 'cannot find 1M.txt -- please contact the developers of pash'
            exit 1
        fi
        "$PASH_TOP/scripts/append_nl_if_not.sh" ./1M.txt
    fi

    if [ ! -f ./10M.txt ]; then
        touch 10M.txt
        for (( i = 0; i < 10; i++ )); do
            cat 1M.txt >> 10M.txt
        done
        "$PASH_TOP/scripts/append_nl_if_not.sh" ./10M.txt
    fi

    if [ ! -f ./100M.txt ]; then
        touch 100M.txt
        for (( i = 0; i < 10; i++ )); do
            cat 10M.txt >> 100M.txt
        done
        "$PASH_TOP/scripts/append_nl_if_not.sh" ./100M.txt
    fi

    if [ ! -f ./1G.txt ]; then
        curl -sf 'http://ndr.md/data/dummy/1G.txt' > 1G.txt
        if [ $? -ne 0 ]; then
            echo 'cannot find 1G.txt -- please contact the developers of pash'
            exit 1
        fi
        "$PASH_TOP/scripts/append_nl_if_not.sh" ./1G.txt
    fi

  # download wamerican-insane dictionary and sort according to machine
  if [ ! -f ./dict.txt ]; then
      curl -sf 'http://ndr.md/data/dummy/dict.txt' | sort > dict.txt
      if [ $? -ne 0 ]; then
          echo 'cannot find dict.txt -- please contact the developers of pash'
          exit 1
      fi
      "$PASH_TOP/scripts/append_nl_if_not.sh" ./dict.txt
    fi

    if [ ! -f ./all_cmds.txt ]; then
        curl -sf 'http://ndr.md/data/dummy/all_cmds.txt' > all_cmds.txt
        if [ $? -ne 0 ]; then
            # This should be OK for tests, no need for abort
            ls /usr/bin/* > all_cmds.txt
        fi
        "$PASH_TOP/scripts/append_nl_if_not.sh" ./all_cmds.txt
    fi


    if [ "$#" -eq 1 ] && [ "$1" = "--full" ]; then
        echo "Generating full-size inputs"
        # FIXME PR: Do we need all of them?

        if [ ! -f ./3G.txt ]; then
            touch 3G.txt
            for (( i = 0; i < 3; i++ )); do
                cat 1G.txt >> 3G.txt
            done
            "$PASH_TOP/scripts/append_nl_if_not.sh" ./3G.txt
        fi

        if [ ! -f ./10G.txt ]; then
            touch 10G.txt
            for (( i = 0; i < 10; i++ )); do
                cat 1G.txt >> 10G.txt
            done
            "$PASH_TOP/scripts/append_nl_if_not.sh" ./10G.txt
        fi

        if [ ! -f ./all_cmdsx100.txt ]; then
            touch all_cmdsx100.txt
            for (( i = 0; i < 100; i++ )); do
                cat all_cmds.txt >> all_cmdsx100.txt
            done
        fi
    fi
}

source_var() {
    if [[ "$1" == "--small" ]]; then
        export IN="$PASH_TOP/evaluation/benchmarks/oneliners/input/small/$2"
        export dict="$PASH_TOP/evaluation/benchmarks/oneliners/input/small/dict.txt"
    else
        export IN="$PASH_TOP/evaluation/benchmarks/oneliners/input/$2"
    fi
}
