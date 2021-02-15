#!/bin/sh


SHELL_TO_JSON_OCAML=/pash/compiler/parser/parse_to_json.native
JSON_TO_SHELL_OCAML=/pash/compiler/parser/json_to_shell.native
JSON_TO_SHELL_C=./json_to_shell2


if [ $# -ne 1 ]
then
    echo "Usage: $0 testFile"
    echo
    exit 1
fi


testFile="$1"


if [ ! -f "$testFile" ]
then
    echo "Error: cannot read '$testFile'!"
    echo
    exit 1
fi


"$SHELL_TO_JSON_OCAML" < "$testFile" > /tmp/json.$$
if [ $? -ne 0 ]
then
    echo "INVALID_INPUT_1: '$testFile' | Unable to run '$SHELL_TO_JSON_OCAML' on '$testFile'"
    exit 1
fi

"$JSON_TO_SHELL_OCAML" < /tmp/json.$$ > /tmp/rt_ocaml.$$
if [ $? -ne 0 ]
then
    echo "INVALID_INPUT_2: '$testFile' | Unable to run '$JSON_TO_SHELL_OCAML' on '/tmp/json.$$'"
    exit 1
fi

"$JSON_TO_SHELL_C" < /tmp/json.$$ > /tmp/rt_c.$$
if [ $? -ne 0 ]
then
    echo "ABORT: '$testFile' | Unable to run '$JSON_TO_SHELL_C' on '/tmp/json.$$'"
    exit 1
fi

diff /tmp/rt_ocaml.$$ /tmp/rt_c.$$
if [ $? -ne 0 ]
then
    diff -w /tmp/rt_ocaml.$$ /tmp/rt_c.$$
    if [ $? -ne 0 ]
    then
        echo "FAIL: '$testFile' | /tmp/json.$$ /tmp/rt_ocaml.$$ /tmp/rt_c.$$"
    else
        echo "FAIL_WHITESPACE: '$testFile' | /tmp/json.$$ /tmp/rt_ocaml.$$ /tmp/rt_c.$$"
    fi
    exit 1
fi

echo "PASS: '$testFile' | /tmp/json.$$ /tmp/rt_ocaml.$$ /tmp/rt_c.$$"
