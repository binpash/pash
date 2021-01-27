#!/bin/sh


SHELL_TO_JSON_OCAML=/pash/compiler/parser/parse_to_json.native
JSON_TO_SHELL_OCAML=/pash/compiler/parser/json_to_shell.native

RT_PY="rt.py"


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

# python3 "$RT_PY" < "$testFile" > /tmp/rt_py.$$
python3 "$RT_PY" "$testFile" > /tmp/rt_py.$$
if [ $? -ne 0 ]
then
    echo "ABORT: '$testFile'"
    exit 1
fi

diff /tmp/rt_ocaml.$$ /tmp/rt_py.$$ > /dev/null
if [ $? -ne 0 ]
then
    diff -w /tmp/rt_ocaml.$$ /tmp/rt_py.$$ > /dev/null
    if [ $? -ne 0 ]
    then
        diff -w /tmp/rt_ocaml.$$ /tmp/rt_py.$$
        echo "FAIL: '$testFile' | /tmp/rt_ocaml.$$ /tmp/rt_py.$$"
    else
        diff /tmp/rt_ocaml.$$ /tmp/rt_py.$$
        echo "FAIL_WHITESPACE: '$testFile' | /tmp/rt_ocaml.$$ /tmp/rt_py.$$"
    fi
    exit 1
fi

echo "PASS: '$testFile' | /tmp/rt_ocaml.$$ /tmp/rt_py.$$"
