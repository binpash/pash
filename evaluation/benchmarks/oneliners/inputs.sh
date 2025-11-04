#!/bin/bash

cd "$(dirname "$0")" || exit 1

#set -e

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

if [ ! -f inputs/1M.txt ]; then


    curl -sfL 'http://atlas-group.cs.brown.edu/data/dummy/1M.txt' >inputs/1M.txt

    if [ $? -ne 0 ]; then
        echo 'cannot find 1M.txt -- please contact the developers of pash'
        exit 1
 
   fi

    "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/1M.txt
fi

if [ ! -f inputs/1G.txt ]; then
    curl -sfL 'http://atlas-group.cs.brown.edu/data/dummy/1G.txt' >inputs/1G.txt

    if [ $? -ne 0 ]; then
        echo 'cannot find 1G.txt -- please contact the developers of pash'
        exit 1
    fi

    "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/1G.txt
fi

# if [ ! -f inputs/3G.txt ]; then

#     touch inputs/3G.txt

#     for i in {1..3}; do
#         cat inputs/1G.txt >>inputs/3G.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/3G.txt
# fi

# if [ ! -f inputs/5G.txt ]; then

#     touch inputs/5G.txt

#     for i in {1..5}; do
#         cat inputs/1G.txt >>inputs/5G.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/5G.txt
# fi

# if [ ! -f inputs/10G.txt ]; then

#     touch inputs/10G.txt

#     for i in {1..10}; do
#         cat inputs/1G.txt >>inputs/10G.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/10G.txt
# fi

# # download wamerican-insane dictionary and sort according to machine
# if [ ! -f inputs/dict.txt ]; then
#     curl -sfL 'http://atlas-group.cs.brown.edu/data/dummy/dict.txt' | sort >inputs/dict.txt

#     if [ $? -ne 0 ]; then
#         echo 'cannot find dict.txt -- please contact the developers of pash'
#         exit 1
#     fi

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/dict.txt
# fi

# if [ ! -f inputs/all_cmds.txt ]; then
#     curl -sfL 'http://atlas-group.cs.brown.edu/data/dummy/all_cmds.txt' >inputs/all_cmds.txt

#     if [ $? -ne 0 ]; then
#         # This should be OK for tests, no need for abort
#         ls /usr/bin/* >inputs/all_cmds.txt
#     fi

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/all_cmds.txt
# fi

