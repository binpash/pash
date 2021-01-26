#!/bin/sh


SHELL_TO_JSON_OCAML=/pash/compiler/parser/parse_to_json.native

PRETTYPRINT_JSON=/pash/compiler/parser/ceda/prettyprint_json

SHELL_TO_JSON_C=/pash/compiler/parser/ceda/parse_to_json2


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


json_ocaml="/tmp/json_ocaml.$$"
json_ocaml_pretty="/tmp/json_ocaml_pretty.$$"
json_c="/tmp/json_c.$$"


"$SHELL_TO_JSON_OCAML" < "$testFile" > "${json_ocaml}"
if [ $? -ne 0 ]
then
    echo "REF_ABORT_1: '$testFile' | Unable to run '$SHELL_TO_JSON_OCAML' on '$testFile'"
    exit 1
fi

"$PRETTYPRINT_JSON" < "${json_ocaml}" > "${json_ocaml_pretty}"
if [ $? -ne 0 ]
then
    echo "REF_ABORT_2: '$testFile' | Unable to run '$PRETTYPRINT_JSON' on '${json_ocaml}'"
    exit 1
fi

"$SHELL_TO_JSON_C" < "$testFile" > "${json_c}"
if [ $? -ne 0 ]
then
    echo "ABORT: '$testFile' | Unable to run '$SHELL_TO_JSON_C' on '$testFile'"
    exit 1
fi

diff "${json_ocaml_pretty}" "${json_c}" > /dev/null
if [ $? -ne 0 ]
then
    diff -w "${json_ocaml_pretty}" "${json_c}" > /dev/null
    if [ $? -ne 0 ]
    then
        diff -w "${json_ocaml_pretty}" "${json_c}"
        echo "FAIL: '$testFile' | ${json_ocaml} ${json_ocaml_pretty} ${json_c}"
    else
        diff "${json_ocaml_pretty}" "${json_c}"
        echo "FAIL_WHITESPACE: '$testFile' | ${json_ocaml} ${json_ocaml_pretty} ${json_c}"
    fi
    exit 1
fi

echo "PASS: '$testFile' | ${json_ocaml} ${json_ocaml_pretty} ${json_c}"
