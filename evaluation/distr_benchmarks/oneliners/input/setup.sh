#!/bin/bash
#set -e

PASH_TOP=${PASH_TOP:-$DISH_TOP/pash}
. "$PASH_TOP/scripts/utils.sh"


# another solution for capturing HTTP status code
# https://superuser.com/a/590170
input_files=("1M.txt" "10M.txt" "100M.txt" "1G.txt" "all_cmds.txt" "all_cmdsx100.txt")
local_fils=("dict.txt")

if [[ "$1" == "-c" ]]; then
    rm -f $input_files "3G.txt" "10G.txt"
    exit
fi

hdfs dfs -mkdir -p /oneliners

if [ ! -f ./1M.txt ]; then
    curl -sf --connect-timeout 10 'atlas-group.cs.brown.edu/data/dummy/1M.txt' > 1M.txt
    if [ $? -ne 0 ]; then
        curl -f 'https://zenodo.org/record/7650885/files/1M.txt' > 1M.txt
        [ $? -ne 0 ] && eexit 'cannot find 1M.txt'
    fi
    append_nl_if_not ./1M.txt
fi

if [ ! -f ./10M.txt ]; then
    touch 10M.txt
    for (( i = 0; i < 10; i++ )); do
        cat 1M.txt >> 10M.txt
    done
fi

if [ ! -f ./100M.txt ]; then
    touch 100M.txt
    for (( i = 0; i < 10; i++ )); do
        cat 10M.txt >> 100M.txt
    done
fi

if [ ! -f ./1G.txt ]; then
    curl -sf --connect-timeout 10 'atlas-group.cs.brown.edu/data/dummy/1G.txt' > 1G.txt
    if [ $? -ne 0 ]; then
        touch 1G.txt
        for (( i = 0; i < 10; i++ )); do
            cat 100M.txt >> 1G.txt
        done
    fi
fi

if [ ! -f ./words ]; then
  curl -sf --connect-timeout 10 'http://ndr.md/data/dummy/words' > words
  if [ $? -ne 0 ]; then
    curl -f 'https://zenodo.org/record/7650885/files/words' > words
    if [ $? -ne 0 ]; then
      if [ $(uname) = 'Darwin' ]; then
        cp /usr/share/dict/web2 words || eexit "cannot find dict file"
      else
        # apt install wamerican-insane
        cp /usr/share/dict/words words || eexit "cannot find dict file"
      fi
    fi
  fi
  append_nl_if_not words
fi

# download wamerican-insane dictionary and sort according to machine
if [ ! -f ./dict.txt ]; then
    curl -sf --connect-timeout 10 'atlas-group.cs.brown.edu/data/dummy/dict.txt' | sort > dict.txt
    if [ $? -ne 0 ]; then
        sort words > sorted_words
    fi
fi

if [ ! -f ./all_cmds.txt ]; then
    curl -sf --connect-timeout 10 'atlas-group.cs.brown.edu/data/dummy/all_cmds.txt' > all_cmds.txt
    if [ $? -ne 0 ]; then
        # This should be OK for tests, no need for abort
        ls /usr/bin/* > all_cmds.txt
    fi
    append_nl_if_not ./all_cmds.txt
fi

if [ ! -f ./all_cmdsx100.txt ]; then
    touch all_cmdsx100.txt
    for (( i = 0; i < 100; i++ )); do
        cat all_cmds.txt >> all_cmdsx100.txt
    done
fi

if [ ! -f ./3G.txt ]; then
    touch 3G.txt
    for (( i = 0; i < 3; i++ )); do
        cat 1G.txt >> 3G.txt
    done
fi
input_files+=("3G.txt")

# Add files with different replication factors
for file in "${input_files[@]}"; do
    hdfs dfs -put $file /oneliners/$file
    rm -f $file
done