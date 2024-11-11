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

# if [ ! -f inputs/10M.txt ]; then
#     touch inputs/10M.txt

#     for (( i = 0; i < 10; i++ )); do
#         cat inputs/1M.txt >>inputs/10M.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/10M.txt
# fi

# if [ ! -f inputs/50M.txt ]; then
#     touch inputs/50M.txt

#     for (( i = 0; i < 50; i++ )); do
#         cat inputs/1M.txt >>inputs/50M.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/50M.txt
# fi

# if [ ! -f inputs/100M.txt ]; then
#     touch inputs/100M.txt

#     for (( i = 0; i < 100; i++ )); do
#         cat inputs/1M.txt >>inputs/100M.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/100M.txt
# fi

# if [ ! -f inputs/200M.txt ]; then
#     touch inputs/200M.txt

#     for (( i = 0; i < 200; i++ )); do
#         cat inputs/1M.txt >>inputs/200M.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/200M.txt
# fi

# if [ ! -f inputs/500M.txt ]; then
#     touch inputs/500M.txt

#     for (( i = 0; i < 500; i++ )); do
#         cat inputs/1M.txt >>inputs/500M.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/500M.txt
# fi

# if [ ! -f inputs/1G.txt ]; then
#     touch inputs/1G.txt

#     for (( i = 0; i < 1000; i++ )); do
#         cat inputs/1M.txt >>inputs/1G.txt
#     done

#     "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/1G.txt
# fi


if [ ! -f inputs/1G.txt ]; then
    curl -sfL 'http://atlas-group.cs.brown.edu/data/dummy/1G.txt' >inputs/1G.txt

    if [ $? -ne 0 ]; then
        echo 'cannot find 1G.txt -- please contact the developers of pash'
        exit 1
    fi

    "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/1G.txt
fi

# download wamerican-insane dictionary and sort according to machine
if [ ! -f inputs/dict.txt ]; then
    curl -sfL 'http://atlas-group.cs.brown.edu/data/dummy/dict.txt' | sort >inputs/dict.txt

    if [ $? -ne 0 ]; then
        echo 'cannot find dict.txt -- please contact the developers of pash'
        exit 1
    fi

    "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/dict.txt
fi

if [ ! -f inputs/3G.txt ]; then

    touch inputs/3G.txt

    for i in {1..3}; do
        cat inputs/1G.txt >>inputs/3G.txt
    done

    "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/3G.txt
fi

if [ ! -f inputs/all_cmds.txt ]; then
    curl -sfL 'http://atlas-group.cs.brown.edu/data/dummy/all_cmds.txt' >inputs/all_cmds.txt

    if [ $? -ne 0 ]; then
        # This should be OK for tests, no need for abort
        ls /usr/bin/* >inputs/all_cmds.txt
    fi

    "$PASH_TOP/scripts/append_nl_if_not.sh" inputs/all_cmds.txt
fi
