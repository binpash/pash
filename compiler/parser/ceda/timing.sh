#!/bin/sh


input_script='/pash/compiler/parser/libdash/ltmain.sh'


if [ $# -eq 1 ]
then
    input_script="$1"
fi


echo "Input script: $input_script"
echo

echo "OCaml:"
time (../parse_to_json.native "$input_script" > /tmp/json.$$; cat /tmp/json.$$ | ../json_to_shell.native | md5sum)
echo

echo "C:"
time (./parse_to_json2 "$input_script" > /tmp/json.$$; cat /tmp/json.$$ | ./json_to_shell2 | md5sum)
echo

echo "Python:"
time (python3 ceda_rt.py "$input_script" | md5sum)
echo

