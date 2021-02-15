#!/bin/sh


SHELL_TO_JSON_OCAML=../parse_to_json.native

PRETTYPRINT_JSON=./prettyprint_json

SHELL_TO_JSON_C=./parse_to_json2


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
    echo "INVALID_INPUT: '$testFile' | Unable to run '$SHELL_TO_JSON_OCAML' on '$testFile'"
    exit 1
fi

"$SHELL_TO_JSON_C" < "$testFile" > "${json_c}"
if [ $? -ne 0 ]
then
    echo "ABORT: '$testFile' | Unable to run '$SHELL_TO_JSON_C' on '$testFile'"
    exit 1
fi


diff "${json_ocaml}" "${json_c}" > /dev/null
if [ $? -ne 0 ]
then
    for f in "${json_ocaml}" "${json_c}"
    do
        "$PRETTYPRINT_JSON" < "${f}" > "${f}.pretty"
        if [ $? -ne 0 ]
        then
            echo "PRETTYPRINT_FAIL: '$testFile' | Unable to run '$PRETTYPRINT_JSON' on '${f}'"
            exit 1
        fi
    done

    diff -w "${json_ocaml}.pretty" "${json_c}.pretty" > /dev/null
    if [ $? -ne 0 ]
    then
        diff -w "${json_ocaml}.pretty" "${json_c}.pretty"
        echo "FAIL: '$testFile' | ${json_ocaml} ${json_c} ${json_ocaml}.pretty ${json_c}.pretty"
    else
        diff "${json_ocaml}" "${json_c}"
        echo "FAIL_WHITESPACE: '$testFile' | ${json_ocaml} ${json_c} ${json_ocaml}.pretty ${json_c}.pretty"
    fi
    exit 1
fi

echo "PASS: '$testFile' | ${json_ocaml} ${json_c} ${json_ocaml}.pretty ${json_c}.pretty"
