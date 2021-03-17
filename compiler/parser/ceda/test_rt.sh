#!/bin/sh


SHELL_TO_JSON_OCAML=../parse_to_json.native
JSON_TO_SHELL_OCAML=../json_to_shell.native

SHELL_TO_JSON_C=./parse_to_json2
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


"$SHELL_TO_JSON_OCAML" < "$testFile" > /tmp/json_ocaml.$$
if [ $? -ne 0 ]
then
    echo "REF_ABORT_1: '$testFile'"
    exit 1
fi

"$JSON_TO_SHELL_OCAML" < /tmp/json_ocaml.$$ > /tmp/rt_ocaml.$$
if [ $? -ne 0 ]
then
    echo "REF_ABORT_2: '$testFile' | /tmp/json_ocaml.$$"
    exit 1
fi

"$SHELL_TO_JSON_C" < "$testFile" > /tmp/json_c.$$
if [ $? -ne 0 ]
then
    echo "ABORT_1: '$testFile'"
    exit 1
fi

"$JSON_TO_SHELL_C" < /tmp/json_c.$$ > /tmp/rt_c.$$
if [ $? -ne 0 ]
then
    echo "ABORT_2: '$testFile' | /tmp/json_c.$$"
    exit 1
fi

diff /tmp/rt_ocaml.$$ /tmp/rt_c.$$
if [ $? -ne 0 ]
then
    diff -w /tmp/rt_ocaml.$$ /tmp/rt_c.$$
    if [ $? -ne 0 ]
    then
        echo "FAIL: '$testFile' | /tmp/rt_ocaml.$$ /tmp/rt_c.$$"
    else
        echo "FAIL_WHITESPACE: '$testFile' | /tmp/rt_ocaml.$$ /tmp/rt_c.$$"
    fi
    exit 1
fi

echo "PASS: '$testFile' | /tmp/rt_ocaml.$$ /tmp/rt_c.$$"
